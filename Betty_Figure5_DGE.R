#reading count tables
library(reshape2)
library(plyr)
library(dplyr)
##read in all tables
tbl<- (1:58278)
#used 58278 bc I know there are 58278 features from featuresCount
for (x in list.files(pattern="*.txt")) {
  u<-read.table(x, header=T,stringsAsFactors = F)
  tbl <- cbind(tbl,u)
}

#get rid of tbl row numbers after if you want
tbl <- tbl[ ,-1]

#truncate column names, or switch to Betty's names
names <- t(read.table('Betty_sampleNames.txt', header = F, stringsAsFactors = F))
colnames(tbl) <- names
tbl$Gene_id<-row.names(tbl)

tbl_mod <- read.table("BettyGeneMatrix.txt", header = F, stringsAsFactors = F)
setwd('/Users/erica/Desktop/Betty_DGEs/')
colnames(tbl_mod)<-tbl_mod[2,]
row.names(tbl_mod) <- tbl_mod[ ,1]
Betty_data_wt <-tbl_mod[4:58281,c(16:29,32:33)]
#convert the matrix data from character to integer values, need to for umap input
Betty_data_wt_N<- type.convert(Betty_data_wt,na.strings = "NA", as.is = FALSE, dec = ".")

#extract just non-invading
Betty_data_wt_NI <- Betty_data_wt_N[ , grepl( "N" , names( Betty_data_wt_N ) ) ] 
#extract just invading
Betty_data_wt_I <- Betty_data_wt_N[ , grepl( "I" , names( Betty_data_wt_N ) ) ] 
#remove controls for now
Betty_data_wt_noC <- Betty_data_wt_N[ ,-(15:16)]

#get list of Gene Ids to convert into gene names
library("biomaRt")
#convertign gene IDs using biomart
Betty_gene_IDs<-rownames(Betty_data_wt_noC)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
E_Ids=Betty_gene_IDs#list of entrez gene ids
Betty_Gene_names <- getBM(attributes=c('ensembl_gene_id_version', 'ensembl_gene_id','external_gene_name'), 
                          filters = 'ensembl_gene_id_version', 
                          values = E_Ids, 
                          mart = ensembl)
annot <- Betty_Gene_names[,c("ensembl_gene_id_version","external_gene_name")]

#remove N1,N2,N4 (due to how they cluster on tsne and umap)
Betty_data_wt_NI_red <- Betty_data_wt_NI[ ,c(3,5,6,7)]
Betty_data_wt_I_red <- Betty_data_wt_I[ ,c(3,5,6,7)]
Betty_data_wt_noC_red <- Betty_data_wt_N[ ,-c(1:4,7,8,15:16,17)]

library(edgeR)
#calculate normalization factors using TMM
Betty_NI_norm = calcNormFactors(Betty_data_wt_NI_red, method='TMM')
Betty_I_norm = calcNormFactors(Betty_data_wt_I_red, method='TMM')
Betty_norm = calcNormFactors(Betty_data_wt_noC_red, method='TMM')
# calculate dispersion factors using:
Betty_NI_disp<- estimateCommonDisp(Betty_data_wt_NI_red, group=NULL, lib.size=NULL, tol=1e-06,
                                   rowsum.filter=5, verbose=FALSE)
Betty_I_disp<- estimateCommonDisp(Betty_data_wt_I_red, group=NULL, lib.size=NULL, tol=1e-06,
                                  rowsum.filter=5, verbose=FALSE)

Betty_disp<- estimateCommonDisp(Betty_data_wt_noC_red, group=NULL, lib.size=NULL, tol=1e-06,
                                rowsum.filter=5, verbose=FALSE)

d <- DGEList(counts=Betty_data_wt_noC_red, group=c(1,2,1,2,1,2,1,2), lib.size=Betty_norm)
Betty_norm = calcNormFactors(Betty_data_wt_noC_red, method='TMM')
Betty_disp<- estimateCommonDisp(Betty_data_wt_noC_red, group=c(1,2,1,2,1,2,1,2), lib.size=Betty_norm, tol=1e-06,
                                rowsum.filter=5, verbose=FALSE)
dim(d) #58278

#filter out lowly expressed genes
keep <- rowSums(cpm(d)>1) >= 0.1
d <- d[keep, , keep.lib.sizes=FALSE]
#test with Fisher exact test
de <- exactTest(d, dispersion=Betty_disp)
de_et_topTags <- topTags(de, n=58278)
#fisher exact test table:
write.table(de_et_topTags, file="ALL_TopTags_et_wt_red.txt", row.names=TRUE, col.names=TRUE, sep="\t")

#testing using GLM
design <- model.matrix(~factor(c(1,2,1,2,1,2,1,2)))
y <- estimateDisp(d,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
qlf_topTags <- topTags(qlf, n=58278)

library(reshape2)
library(plyr)
library(dplyr)
library(reshape)

#get individualcpm normalized values for each sample to make a heatmap
qlf_topTags_ind <-rownames(topTags(qlf, n=58278)$table)
tagsbysample_qlf<-as.data.frame(cpm(d)[qlf_topTags_ind, order(d$samples$group)])
tagsbysample_qlf$Gene_id <-rownames(tagsbysample_qlf)

#get gene names
annot <- read.table("ALL_TopTags_qlf_wt_red_GeneNames.txt", header = T, stringsAsFactors = F)

#retrieve genes with a significant pvalue
qlf_topTags_id<-subset(qlf_topTags$table, FDR < 0.05) #244
qlf_topTags_id$Gene_id <-rownames(qlf_topTags_id)
qlf_heatmap_mat <- merge(qlf_topTags_id,annot, by.x="Gene_id",by.y="GenestableIDversion") #now 216
qlf_heatmap_mat <- qlf_heatmap_mat[ ,c(1,7)]
qlf_heatmap_mat_ind <- merge(qlf_heatmap_mat,tagsbysample_qlf,by="Gene_id")

##facetted heatmap
qlf_heatmap_mat_ind <- qlf_heatmap_mat_ind[ ,-1]
qlf_heatmap_mat_ind_melt <- melt(qlf_heatmap_mat_ind,id.vars="Genename") 
qlf_heatmap_mat_ind_melt_split <- t(data.frame(strsplit(as.character(qlf_heatmap_mat_ind_melt$variable),"",fixed = T)))
qlf_heatmap_mat_ind_mod <- cbind(qlf_heatmap_mat_ind_melt,qlf_heatmap_mat_ind_melt_split)
#qlf_heatmap_mat_ind_mod$gene_name<-qlf_heatmap_mat_ind_mod$external_gene_name
#rownames(qlf_heatmap_mat_ind_mod) <- qlf_heatmap_mat_ind_mod$gene_name
qlf_heatmap_mat_ind_mod_red <- qlf_heatmap_mat_ind_mod[ ,c(1,3,4,5)]
colnames(qlf_heatmap_mat_ind_mod_red) <- c('gene','value','state','exp')

##making a gene clustered, viridis, facetted heat map for qlf signficant genes
library(ggplot2)
library(viridis)
library(ggExtra)
col_breaks <- c(1,5,10,20,50,100,250,500,750,1000)

##cluster data manually to get order on the geom_tile()
rownames(qlf_heatmap_mat_ind) <- qlf_heatmap_mat_ind$Genename
qlf_heatmap_mat_ind <- qlf_heatmap_mat_ind[ ,-1]
hr_qlf <- hclust(as.dist(1-cor(t(qlf_heatmap_mat_ind), method="pearson")),
                 method="average")

gene_order <- hr_qlf$labels

pdf(file="Betty_facet_qlf_heat_mod.pdf",width=4, height=10)
par(mar=c(7,4,4,2)+0.1) 
ggplot(qlf_heatmap_mat_ind_mod_red,aes(state,gene,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis_c(name="Gene \nExpression (CPM)",limits = c(0,50), breaks =col_breaks,option="viridis",
                       na.value = "yellow") + 
  scale_y_discrete(limits=gene_order, labels=gene_order) +
  scale_x_discrete(limits=c("N","I"), labels=c("N", "NI")) +
  facet_grid(.~exp, switch = "x", scales = "free_x", space = "free_x") +
  labs(x="State", y="Gene") +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 16, face="bold"))+
  theme(axis.text.y=element_text(size=5)) +
  theme(strip.background = element_rect(colour="white"))+
  theme(plot.title=element_text(hjust=0))+
  theme(axis.ticks=element_blank())+
  theme(axis.text.y=element_text(size=3))+
  theme(axis.text.x=element_text(size=12))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=6))+
  theme(strip.text.x = element_text(size=0))
  guides(fill=guide_colorbar(nbin=20)) +
  removeGrid() 
dev.off()

hmap_order <- data.frame(qlf_heatmap_mat_ind[(hr_qlf$labels), ])
write.csv(hmap_order ,"Betty_qlfheatmap.csv")






