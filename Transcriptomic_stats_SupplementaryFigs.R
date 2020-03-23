##Reads vs gene counts, S10A##
#retrieve gene table
GeneMat <- read.table("BettyGeneMatrix.txt", header = F, stringsAsFactors = F)
#subset wildtype and controls
colnames(GeneMat) <- GeneMat[2,]
rownames(GeneMat)<-GeneMat$Sample
Betty_data_wt <-GeneMat[4:58281,c(1,16:29,32:33)]

#get gene lengths
gene_l <- read.table("Human_gene_length.txt", header = T, stringsAsFactors = F)
#getting gene lengths from chr end-chr start
gene_l$gene_L <- (gene_l$Gene_end-gene_l$Gene_start)
#This file^ actually has gene start,emd, transcript length, GC content, gene type information
#retrieved from Biomart
#getting only what I need for right now, Gene lengths:
gene_l_red <- gene_l[ ,c("Gene_name", "GeneStableIDversion","gene_L")]
gene_l_red_u <- unique(gene_l_red)
#Merge with gene lengths (Merging anf preserving order to look for duplicates and preserve order)
Betty_data_wt$id <- 1:nrow(Betty_data_wt) 
Betty_data_wt_L <- merge(Betty_data_wt,gene_l_red_u, by.x="Sample", by.y="GeneStableIDversion", sort=F)
ordered <- Betty_data_wt_L[order(Betty_data_wt_L$id), ]
#removing id and GeneStableIDversion columns
Betty_data_wt_preMelt <-ordered[ ,-c(1,18)] 
#melt this matrix
Betty_data_wt_melt <- melt(Betty_data_wt_preMelt,id.vars= c("Gene_name","gene_L"))
Betty_data_wt_melt_N <- cbind(Betty_data_wt_melt[ ,1:3],type.convert(Betty_data_wt_melt$value,na.strings = "NA", as.is = FALSE, dec = "."))  
colnames(Betty_data_wt_melt_N) <- c("Gene_name","gene_L","variable","cpm")
#remove genes with CPM less than 0.1
Betty_data_wt_melt_robust<-subset(Betty_data_wt_melt_N,cpm > 0.1)

#get cumulative CPM and gene length data for each sample
Betty_data_wt_melt_robust_sum<- ddply(Betty_data_wt_melt_robust, c("variable"), summarise,
                                      sum = sum(cpm), gene_L=sum(gene_L),gene_N=length(Gene_name))

##merge this with read counts
#Get number of reads per sample
reads_sample <- read.table("Sample_reads.txt", header = T, stringsAsFactors = F)
Betty_data_stats <- merge(Betty_data_wt_melt_robust_sum,reads_sample,by.x="variable", by.y="sample_name")
Betty_data_sample.names <- Betty_data_stats$variable

Betty_data_stats$cond[grepl( "N" , Betty_data_stats$variable)]<-"N"
Betty_data_stats$cond[grepl( "I" , Betty_data_stats$variable )]<-"I"
Betty_data_stats$cond[grepl( "wpWT" , Betty_data_stats$variable )]<-"Control well"
Betty_data_stats$cond[grepl( "gWT" , Betty_data_stats$variable )]<-"Control gel"

##plot all these:
library(ggplot2)
pdf("Betty_geneVreads_trans.pdf")
ggplot(data=Betty_data_stats) + geom_point(aes(x=number_of_reads_mapped,y=gene_N,colour=cond,size=gene_L)) +
  scale_colour_manual("Condition", values = c("#c11717","#F8766D","#00BA38","#619CFF"),
                      labels=c("Control Well","Control Gel","Invading Cells","Non-invading Cells")) +
  scale_size_continuous("Length of \nTranscriptome (bp) \ncovered",range = c(5,8)) +
  ylab("Number of Genes Detected") + xlab("Number of Reads") +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title.align = 0.5) +
  theme(axis.text.y = element_text(colour="black", size = 18)) +
  theme(axis.text.x = element_text(colour="black", size = 18, angle=25, hjust=1)) +
  theme(axis.title = element_text(colour="black", size = 20, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=8))) 
dev.off()

##UMAP_S10B##
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
setwd('/Users/erica/Desktop/Betty_DGE/Betty_stats/')
names <- t(read.table('Betty_sampleNames.txt', header = F, stringsAsFactors = F))
setwd('/Users/erica/Desktop/Betty_DGE/')
colnames(tbl) <- names
tbl$Gene_id<-row.names(tbl)

tbl <- tbl[ ,-33]
Betty_data_t <- t(tbl)
Betty_data_sample.names <- rownames(Betty_data_t) 
Betty_data_gene.names <- colnames(Betty_data_t)

library(umap)
Betty_umap_t <- umap(Betty_data_t,n_epochs=500)
head(Betty_umap_t$layout, 3)
Betty_umap_t_lay <- as.data.frame(Betty_umap_t$layout) 
colnames(Betty_umap_t_lay)<- c('Dim1','Dim2')
Betty_umap_t_lay$sample <- rownames(Betty_umap_t_lay)

#umap plot with labels and colours
pdf(file="Betty_umap.pdf",width=7, height=4)
ggplot(Betty_umap_t_lay, aes(Dim1,Dim2,colour=sample)) +
  geom_point() +
  scale_fill_manual(values = my.cols,
                    labels=Betty_data_sample.names) +
  geom_text(aes(label=Betty_data_sample.names),hjust=1, vjust=1)
dev.off()

##Control versus CIMMS heatmap, S10C##
GeneMat<- read.table("BettyGeneMatrix.txt", header = F, stringsAsFactors = F)
#subset wildtype and controls
colnames(GeneMat) <- GeneMat[2,]
rownames(GeneMat)<-GeneMat$Sample
Betty_data_wt <-GeneMat[4:58281,c(1,16:29,32:33)]
#add Gene names
gene_l <- read.table("Human_Gene_length.txt", header = T, stringsAsFactors = F)
gene_names <- gene_l[ ,c("Gene_name", "GeneStableIDversion")]
gene_names_u <- unique(gene_names)
#Merge tp attach gene names
Betty_data_wt_names <- merge(Betty_data_wt,gene_names_u, by.x="Sample", by.y="GeneStableIDversion", sort=F)
##melt and split to parse out labels, prep for facetted heatmap
Betty_data_wt_names_mat <- Betty_data_wt_names[ ,-1]
row.names(Betty_data_wt_names_mat) <- Betty_data_wt_names[ ,1]
#Betty_data_wt_names_mat <- type.convert(Betty_data_wt_names_mat,na.strings = "NA", as.is = FALSE, dec = ".")
Betty_data_wt_names_mat_melt <- melt(Betty_data_wt_names_mat,id.vars=c("Gene_name")) 
Betty_data_wt_names_mat_melt_split_list <- (strsplit(as.character(Betty_data_wt_names_mat_melt$variable),"",fixed = T))
##had to make this into a list bc keep getting error of unequal rows 2,3,4...
Betty_data_wt_names_mat_melt_split <- t(data.frame(lapply(Betty_data_wt_names_mat_melt_split_list, `length<-`, max(lengths(Betty_data_wt_names_mat_melt_split_list)))))
Betty_data_wt_names_mat_mod <- cbind(Betty_data_wt_names_mat_melt,Betty_data_wt_names_mat_melt_split)
Betty_data_wt_names_mat_mod_red <- Betty_data_wt_names_mat_mod[ ,c(1,3,4,5)]
colnames(Betty_data_wt_names_mat_mod_red) <- c('gene','value','state','exp')
Betty_data_wt_names_mat[ ,1:16] <- type.convert(Betty_data_wt_names_mat[ ,1:16],na.strings = "NA", as.is = FALSE, dec = ".")
##sum the N and Is from each experiment to simulate a control
Betty_data_wt_names_mat$Exp1<-(Betty_data_wt_names_mat$N1 + Betty_data_wt_names_mat$I1)
Betty_data_wt_names_mat$Exp2<-(Betty_data_wt_names_mat$N2 + Betty_data_wt_names_mat$I2)
Betty_data_wt_names_mat$Exp3<-(Betty_data_wt_names_mat$N3 + Betty_data_wt_names_mat$I3)
Betty_data_wt_names_mat$Exp4<-(Betty_data_wt_names_mat$N4 + Betty_data_wt_names_mat$I4)
Betty_data_wt_names_mat$Exp5<-(Betty_data_wt_names_mat$N5 + Betty_data_wt_names_mat$I5)
Betty_data_wt_names_mat$Exp6<-(Betty_data_wt_names_mat$N6 + Betty_data_wt_names_mat$I6)
Betty_data_wt_names_mat$Exp7<-(Betty_data_wt_names_mat$N7 + Betty_data_wt_names_mat$I7)
#make a matrix with only sums and controls
Betty_data_wt_mat_calc <- Betty_data_wt_names_mat[ ,c(15,16,18:24)]

##make heatmap
#cluster rows in matrix
qlf_mat <- qlf_heatmap[ ,-1]
row.names(qlf_mat)<-qlf_heatmap[ ,1] 
hr_control <- hclust(as.dist(1-cor(t(qlf_mat), method="pearson")),
                     method="average")
gene_order <- hr_control$labels
#cluster samples for fun, but will not use for this heatmap
hc_control <- hclust(as.dist(1-cor(qlf_mat, method="spearman")), method="average")

qlf_mat<- type.convert(qlf_mat,na.strings = "NA", as.is = FALSE, dec = ".")

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("Blue", "white", "Red"))(n = 14)
#n must always be 1 less than the length(col_breaks)
#making the heatmap
library(gplots)
library(RColorBrewer)
library(svDialogs)
library(viridis)
qlf_mat <- as.matrix(qlf_mat)
qlf_mat_log <- log(qlf_mat)

options("scipen"=100, "digits"=4)
lmat = rbind(c(0,4),c(0,3),c(2,1))
lwid = c(0.5,4)
lhei = c(1,0.5,4)

pdf(file='Betty_ctrl_heatmap.pdf', width=5, height=8)
par(mar=c(7,4,4,2)+0.1) 
#you can play with col_breaks and colorRampPalette above to change colour saturation of the heatmap
#col_breaks <- c(50,100,200,300,500,1000,2000,5000,7500,10000,15000,30000,50000,100000,10000000)
col_breaks <- c(0,1,2,3,5,7,10,12,14,16,18,20,25,30)
heatmap.2(qlf_mat_log,    # data matrix
          trace="none", 
          #labRow = NULL,
          margins =c(5,1),     # widens margins around plot
          col=viridis,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          Colv=as.dendrogram(hc_control),
          Rowv=as.dendrogram(hr_control),
          hclustfun = hclust,
          keysize = 1,
          labRow = FALSE,
          key.xlab = "log(cpm)",
          key.title = NA,
          densadj = 0.5,
          density.info="density", 
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)                # turn off column clustering
dev.off()  # close the PNG device

hmap_order <- data.frame(qlf_mat[rev(hr_control$labels[hr_control$order]), hc_control$labels[hc_control$order]])
write.csv(hmap_order ,"Betty_ctrl_qlf_heatmap.csv")
