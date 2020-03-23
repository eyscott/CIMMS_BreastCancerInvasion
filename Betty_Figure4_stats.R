##packages needed for this workflow
library(reshape2)
library(plyr)
library(dplyr)

#plot length of genes vs # of reads
#Get gene matrix
GeneMat <- read.table("BettyGeneMatrix.txt", header = F, stringsAsFactors = F)

#subset wildtype and controls
colnames(GeneMat) <- GeneMat[2,]
rownames(GeneMat)<-GeneMat$Sample
Betty_data_wt <-GeneMat[4:58281,c(1,16:29,32:33)]

#get gene lengths
gene_l <- read.table("Human_gene_length.txt", header = T, stringsAsFactors = F) #retrieved from Biomart
#getting gene lengths from chr end -chr start
gene_l$gene_L <- (gene_l$Gene_end-gene_l$Gene_start)
#This file^ actually has gene start,emd, transcript length, GC content, gene type information
#Retrieving gene lengths:
gene_l_red <- gene_l[ ,c("Gene_name", "GeneStableIDversion","gene_L")]
gene_l_red_u <- unique(gene_l_red)
#Merge with gene lengths (Merging and preserving order to look for duplicates)
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
setwd('/Users/erica/Desktop/Betty_DGEs')
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
  #geom_text(aes(x=number_of_reads_mapped,y=gene_N,label=Betty_data_sample.names),hjust=1, vjust=1,size=8) +
  #theme(plot.margin = margin(0, 0, 1, 2, "cm")) 
dev.off()

######making a heatmap for control vs control vs sum of exps#########
setwd('/Users/erica/Desktop/Betty_DGEs/Betty_stats/')
#reading count tables
library(reshape2)
library(plyr)
library(dplyr)
#make a matrix with sum of exps#
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
Betty_data_wt_mat_calc <- type.convert(Betty_data_wt_mat_calc,na.strings = "NA", as.is = FALSE, dec = ".")
write.table(Betty_data_wt_mat_calc,"CtrlGeneMat.txt")
library(edgeR)
Betty_norm = calcNormFactors(Betty_data_wt_mat_calc, method='TMM')
Betty_disp<- estimateCommonDisp(Betty_data_wt_mat_calc, group=c(1,1,2,2,2,2,2,2,2), lib.size=Betty_norm, tol=1e-06,
                                rowsum.filter=5, verbose=FALSE)
d <- DGEList(counts=Betty_data_wt_mat_calc, group=c(1,1,2,2,2,2,2,2,2), lib.size=Betty_norm)
dim(d) #58278
#testing using GLM
design <- model.matrix(~factor(c(1,1,2,2,2,2,2,2,2)))
y <- estimateDisp(d,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
qlf_topTags <- topTags(qlf, n=58278)
#GLM results:
write.table(qlf_topTags, file="CvsExp_TopTags_qlf.txt", row.names=TRUE, col.names=TRUE, sep="\t")
#retrieve significant genes
qlf_topTags_ind <-rownames(topTags(qlf, n=55067)$table)
tagsbysample_qlf<-as.data.frame(cpm(d)[qlf_topTags_ind, order(d$samples$group)])
tagsbysample_qlf$Gene_id <-rownames(tagsbysample_qlf)
#retrieve genes with a significant pvalue
qlf_topTags_id<-subset(qlf_topTags$table, PValue < 0.01) #441
qlf_topTags_id$Gene_id <-rownames(qlf_topTags_id)
#merge with CPM gene matrix
#qlf_mat <- merge(qlf_heatmap_mat,tagsbysample_qlf,by="Gene_id")

#attach external gene names
qlf_heatmap_mat <- merge(qlf_topTags_id,gene_names_u, by.x="Gene_id", by.y="GeneStableIDversion", sort=F)
qlf_heatmap_mat_ind <- merge(qlf_heatmap_mat,tagsbysample_qlf,by="Gene_id")
qlf_heatmap_mat_CPM_name <- qlf_heatmap_mat_ind[ ,c(7:16)]
write.table(qlf_heatmap_mat_CPM_name, "qlfCtrlGeneMat.txt")
setwd('/Users/erica/Desktop/Betty_DGEs')
qlf_heatmap <- read.table("qlfCtrlGeneMat.txt",header = T, stringsAsFactors = F)
#make the plots
library(ggplot2)
library(viridis)
library(ggExtra)
col_breaks <- c(seq(0,100,10), seq(1000,1000000000,1000))

#cluster rows in matrix
qlf_mat <- qlf_heatmap[ ,-1]
hr_control <- hclust(as.dist(1-cor(t(qlf_mat), method="pearson")),
                     method="average")
gene_order <- hr_control$labels
#cluster samples for fun, but will not use for this heatmap
hc_control <- hclust(as.dist(1-cor(qlf_mat, method="spearman")), method="average")

row.names(qlf_mat)<-qlf_heatmap[ ,1] 
qlf_mat<- type.convert(qlf_mat,na.strings = "NA", as.is = FALSE, dec = ".")

#just need to melt for plotting
qlf_heatmap_mat_CPM_name_melt <- melt(qlf_heatmap_mat_CPM_name,id.vars = "Gene_name")

pdf(file='Betty_Controls_qlf_heat.png', width=10, height=8)
par(mar=c(7,4,4,2)+0.1) 
ggplot(qlf_heatmap_mat_CPM_name_melt,aes(variable,Gene_name,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis_c(name="Gene \nExpression",
                       breaks=col_breaks,option="viridis") + 
  scale_y_discrete(limits=gene_order,labels=gene_order) +
  labs(x="Experiment", y="Gene") +
  theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 14))+
  theme(axis.text.y=element_text(size=6)) +
  theme(strip.background = element_rect(colour="white"))+
  theme(plot.title=element_text(hjust=0))+
  theme(axis.ticks=element_blank())+
  theme(axis.text.y=element_text(size=3))+
  theme(axis.text.x=element_text(size=10))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=6))+
  guides(fill=guide_colorbar(nbin=20)) 
  removeGrid()#ggExtra
dev.off()

##getting a table with the genes in order from above
hmap_order <- data.frame(qlf_heatmap_mat_CPM_name[rev(hr_control$labels[hr_control$order]), hc_control$labels[hc_control$order]])
write.csv(hmap_order ,"Betty_control_heatmap.csv")


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


##umap
##did this on HPC
umap_data <- read.table("Betty_umap.txt", header=T, stringsAsFactors = F)
umap_data$cond[grepl( "N" , umap_data$sample)]<-"N"
umap_data$cond[grepl( "I" , umap_data$sample)]<-"I"
umap_data$cond[grepl( "WT" , umap_data$sample)]<-"Cont"
Betty_data_sample.names <-  umap_data$sample

library(ggplot2)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
pdf(file="Betty_wt_cond.pdf",width=10, height=8)
ggplot(umap_data, aes(Dim1,Dim2,colour=cond)) +
  geom_point() +
  scale_colour_manual("Condition", values = c("#F8766D","#00BA38","#619CFF"),
                      labels=c("Control","Invading Cells","Non-invading Cells")) +
  geom_text(aes(label=Betty_data_sample.names),hjust=1, vjust=1,show.legend = FALSE, size=4) +
  theme(legend.title = element_text(colour="black", size = 18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text.y = element_text(colour="black", size = 18)) +
  theme(axis.text.x = element_text(colour="black", size = 18)) +
  theme(axis.title.y = element_text(colour="black",size = 24)) +
  theme(axis.title.x = element_text(colour="black",size = 24)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_size(range = c(2,12)) +
  guides(colour = guide_legend(override.aes = list(size=8))) 
 dev.off()

