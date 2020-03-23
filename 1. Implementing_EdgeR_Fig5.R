#prerequisites for formatting data
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

##read in all tables
tbl<- (1:58278)
#used 58278 bc I know there are 58278 features from featuresCount
for (x in list.files(pattern="*.txt")) {
  u<-read.table(x, header=T,stringsAsFactors = F)
  tbl <- cbind(tbl,u)
}

#get rid of tbl row numbers after if you want
tbl <- tbl[ ,-1]

#truncate column names by replacing with manual-made column names
names <- t(read.table('Betty_sampleNames.txt', header = F, stringsAsFactors = F))
setwd('/Users/erica/Desktop/Betty_DGEs/')
colnames(tbl) <- names
tbl$Gene_id<-row.names(tbl)

write.table(tbl,"BettyGeneMatrix.txt")
#manually added treatment groups (as rows) in excel, sorry :| :)

tbl_mod <- read.table("BettyGeneMatrix.txt", header = F, stringsAsFactors = F)
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
#converting gene IDs using biomart
Betty_gene_IDs<-rownames(Betty_data_wt_noC)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
E_Ids=Betty_gene_IDs#list of entrez gene ids
Betty_Gene_names <- getBM(attributes=c('ensembl_gene_id_version', 'ensembl_gene_id','external_gene_name'), 
                          filters = 'ensembl_gene_id_version', 
                          values = E_Ids, 
                          mart = ensembl)
annot <- Betty_Gene_names[,c("ensembl_gene_id_version","external_gene_name")]

library(edgeR)
#calculate normalization factors using TMM
d <- DGEList(counts=Betty_data_wt_noC_red, group=c(1,2,1,2,1,2,1,2), lib.size=Betty_norm)
Betty_norm = calcNormFactors(Betty_data_wt_noC_red, method='TMM')
#Betty_norm = 1.1949235 0.9554808 1.0388561 0.9758201 1.1799904 1.0964608 1.0482861 0.6370322
Betty_disp<- estimateCommonDisp(Betty_data_wt_noC_red, group=c(1,2,1,2,1,2,1,2), lib.size=Betty_norm, tol=1e-06,
                                rowsum.filter=5, verbose=FALSE) 
# Betty_disp = 1.2

dim(d) #[1] 58278     8

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

#get individual cpm normalized values for each sample to make a heatmap
qlf_topTags_ind <-rownames(topTags(qlf, n=58278)$table)
tagsbysample_qlf<-as.data.frame(cpm(d)[qlf_topTags_ind, order(d$samples$group)])
tagsbysample_qlf$Gene_id <-rownames(tagsbysample_qlf)

#get gene names
annot2 <- read.table("ALL_TopTags_qlf_wt_red_GeneNames.txt", header = T, stringsAsFactors = F) #got from Biomart

#retrieve genes with a significant FDR (FDR < 0.05)
qlf_topTags_id<-subset(qlf_topTags$table, FDR < 0.05) #244 genes
qlf_topTags_id$Gene_id <-rownames(qlf_topTags_id)
qlf_heatmap_mat <- merge(qlf_topTags_id,annot2, by.x="Gene_id",by.y="GenestableIDversion") #now 216, aka, only 216 have gene names
qlf_heatmap_mat <- qlf_heatmap_mat[ ,c(1,7)]
qlf_heatmap_mat_ind_names <- merge(qlf_heatmap_mat,tagsbysample_qlf,by="Gene_id")

## Making FIGURE 5 heatmap
##facetted heatmap
qlf_heatmap_mat_ind <- qlf_heatmap_mat_ind_names[ ,-1]
qlf_heatmap_mat_ind_melt <- gather(qlf_heatmap_mat_ind,"variable", "value", 2:9)
qlf_heatmap_mat_ind_melt_split <- t(data.frame(strsplit(as.character(qlf_heatmap_mat_ind_melt$variable),"",fixed = T)))
qlf_heatmap_mat_ind_mod <- cbind(qlf_heatmap_mat_ind_melt,qlf_heatmap_mat_ind_melt_split)
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

pdf(file="Betty_separate_exp_facetted.pdf",width=4, height=10)
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

## Make average columns for Invading (I) versus Non-Invading (NI)
#extract just non-invading
Betty_qlf_heatmap_mat_NI <- cbind(Gene_name=qlf_heatmap_mat_ind_names$Genename,qlf_heatmap_mat_ind_names[ ,grepl("N",colnames(qlf_heatmap_mat_ind_names))])
#extract just invading
Betty_qlf_heatmap_mat_I <- cbind(Gene_name=qlf_heatmap_mat_ind_names$Genename,qlf_heatmap_mat_ind_names[ ,grepl("I",colnames(qlf_heatmap_mat_ind_names))]) 

# "gather" to average
Betty_qlf_heatmap_I_melt <- gather(Betty_qlf_heatmap_mat_I,"variable", "value",2:5)
Betty_qlf_heatmap_N_melt <- gather(Betty_qlf_heatmap_mat_NI,"variable", "value", 2:5)
# then employ dpylr for averaging
Betty_qlf_heatmap_I_melt_mean<- ddply(Betty_qlf_heatmap_I_melt, c("Gene_name"), summarise,
                                      mean_I=mean(value))
Betty_qlf_heatmap_N_melt_mean<- ddply(Betty_qlf_heatmap_N_melt, c("Gene_name"), summarise,
                                      mean_N=mean(value))
##attach these
Betty_mean_N_and_I <- cbind(Betty_qlf_heatmap_I_melt_mean,Betty_qlf_heatmap_N_melt_mean)
#now make the logFC column
Betty_mean_FC_I <- (Betty_mean_N_and_I$mean_I-Betty_mean_N_and_I$mean_N)
Betty_mean_N_and_I_FC<- cbind(Betty_mean_N_and_I[ ,c(1,2,4)],Betty_mean_FC_I )
#prep for heatmap, make into matrix
rownames(Betty_mean_N_and_I_FC) <- Betty_mean_N_and_I_FC$Gene_name
Betty_mean_N_and_I_FC_mat <- as.matrix(type.convert(Betty_mean_N_and_I_FC[ ,c(2,3,4)],na.strings = "NA", as.is = FALSE, dec = "."))
#cluster the rows
hr <- hclust(as.dist(1-cor(t(Betty_mean_N_and_I_FC_mat[ ,c(1,2)]), method="pearson")),
             method="average")
Betty_mean_N_and_mat <- Betty_mean_N_and_I_FC_mat[ ,1:2]
setwd('/Users/erica/Desktop/Betty_DGE/Figs')
pdf(file='Betty_F5A.pdf', width=5, height=8,bg="white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
col_breaks <-c(seq(0,2.4,0.2),seq(3,5,0.5),6,7.5,10,12.5,15)
heatmap.2(Betty_mean_N_and_mat,    # data matrix
          trace="none", 
          margins =c(8,8),     # widens margins around plot
          col=viridis,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          labCol = c("Average I","Average N"),
          cexCol = 1.5,
          dendrogram="none",     # only draw a row dendrogram
          Colv=F,
          Rowv=as.dendrogram(hr),
          hclustfun = hclust,
          keysize = 1,
          labRow = FALSE,
          key.xlab = "CPM",
          key.title = NA,
          densadj = 0.5,
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)                # turn off column clustering
dev.off()  # close the PNG device

hc <- hclust(as.dist(1-cor(Betty_mean_N_and_I_FC_mat, method="spearman")), method="average")
hmap_order <- data.frame(Betty_mean_N_and_I_FC_mat[rev(hr$labels[hr$order]), hc$labels[hc$order]])
write.csv(hmap_order ,"Betty_F5A_heatmap.csv")

