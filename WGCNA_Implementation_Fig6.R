#preprocess BettyGeneMatrix with just wt cell line and conventional gene names
WGCNA_in1 <- read.table("Betty_stats/BettyGeneMatrix.txt", header = F, stringsAsFactors = F)
colnames(WGCNA_in1)<-WGCNA_in1[2,]
#make a metadata sheet for later
WGCNA_meta <- WGCNA_in1[1:3,c(16:29,32:33)]
WGCNA_in1 <-WGCNA_in1[4:58281,c(16:29,34)]
##keep only 3,5,6,7 (like DE analysis)
WGCNA_in1.5<- WGCNA_in1[ ,c(5,6,9:15)]

#remove love expressed genes
library(reshape2)
library(plyr)
library(dplyr)
library(reshape)
WGCNA_in1_melt <- melt(WGCNA_in1, id.vars= c("Gene_id"))
WGCNA_in1_melt$value <- as.numeric(as.character(WGCNA_in1_melt$value))
WGCNA_in1_sum<- ddply(WGCNA_in1_melt, c("Gene_id"), summarise,
                sum = sum(value), sd = sd(value),
                sem = sd(value)/sqrt(length(value)))
WGCNA_in1_sum_merge <- merge(WGCNA_in1,WGCNA_in1_sum,by="Gene_id")
WGCNA_in1_high <- subset(WGCNA_in1_sum_merge, c(sum > 10))
#variance stabilize remaining data
row.names(WGCNA_in1_high) <- WGCNA_in1_high$Gene_id
WGCNA_in1_high_N<- type.convert(WGCNA_in1_high[ ,c(2:15)],na.strings = "NA", as.is = FALSE, dec = ".")
WGCNA_in1_high_log <- log(WGCNA_in1_high_N+1)

library(reshape2)
library(plyr)
library(dplyr)
library(reshape)
WGCNA_in1.5_melt <- melt(WGCNA_in1.5, id.vars= c("Gene_id"))
WGCNA_in1.5_melt$value <- as.numeric(as.character(WGCNA_in1.5_melt$value))
WGCNA_in1.5_sum<- ddply(WGCNA_in1.5_melt, c("Gene_id"), summarise,
                      sum = sum(value), sd = sd(value),
                      sem = sd(value)/sqrt(length(value)))
WGCNA_in1.5_sum_merge <- merge(WGCNA_in1.5,WGCNA_in1.5_sum,by="Gene_id")
WGCNA_in1.5_high <- subset(WGCNA_in1.5_sum_merge, c(sum > 10))
#variance stabilize remaining data
row.names(WGCNA_in1.5_high) <- WGCNA_in1.5_high$Gene_id
WGCNA_in1.5_high_N<- type.convert(WGCNA_in1.5_high[ ,c(2:9)],na.strings = "NA", as.is = FALSE, dec = ".")
WGCNA_in1.5_high_log <- log(WGCNA_in1.5_high_N+1)

library(WGCNA)
options(stringsAsFactors = FALSE)
WGCNA_in2 = as.data.frame(t(WGCNA_in1_high_log))
#check for genes and samples with too many missing values
gsg = goodSamplesGenes(WGCNA_in2, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(WGCNA_in2)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(WGCNA_in2)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  WGCNA_in2 = WGCNA_in2[gsg$goodSamples, gsg$goodGenes]
}

#cluster samples to check for outliers
sampleTree = hclust(dist(WGCNA_in2), method = "average");
sizeGrWindow(12,9)
pdf(file = "BettySampleCluster_outliers.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
graphics.off()
#remove outliers by height
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 300, minSize = 2)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==c(1))
WGCNA_in3 = WGCNA_in2[keepSamples, ]
nGenes = ncol(WGCNA_in3)
nSamples = nrow(WGCNA_in3)
##now bring in the trait data
WGCNA_meta_t <- as.data.frame(t(WGCNA_meta))
colnames(WGCNA_meta_t)<-c('trt','exp','num')
WGCNA_meta_t$wgcna[grepl( "I" , WGCNA_meta_t$exp)]<- 1
WGCNA_meta_t$wgcna[grepl( "N" , WGCNA_meta_t$exp)]<-2
WGCNA_meta_t$wgcna[grepl( "g" , WGCNA_meta_t$exp)]<-3
WGCNA_meta_t$wgcna[grepl( "wp" , WGCNA_meta_t$exp)]<-4
#manually remove samples removed bc of outliers I3,N2,N6
#manually remove controls
WGCNA_meta_2 <- WGCNA_meta_t[-c(15,16),]
# Re-cluster samples
sampleTree2 = hclust(dist(WGCNA_in3), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(WGCNA_meta_2$wgcna, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    #groupLabels = rownames(WGCNA_meta_t),
                    main = "Sample dendrogram and trait heatmap")
save(WGCNA_in3, WGCNA_meta_2, file = "Betty_wgcna_dataInput.RData")

allowWGCNAThreads()
lnames = load(file = "Betty_wgcna_dataInput.RData")
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_in3, powerVector = powers, verbose = 10)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

save(WGCNA_in3, WGCNA_meta_2, sft, file = "2_Betty_wgcna_postPower.RData")

#####do this on HPC
lnames = load(file = "2_Betty_wgcna_postPower.RData")
lnames
softPower = 10;
adjacency = adjacency(WGCNA_in3, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
save(TOM, dissTOM, file = "Tom_Betty.RData")
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(WGCNA_in3, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.08
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(WGCNA_in3, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
pdf(file = "geneDendro_3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Betty_networkConstruction_stepByStep.RData")

##########################end of HPC usage


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Betty_wgcna_dataInput.RData")
lnames = load(file = "Betty_ME.RData")
#####
MEs <- MEs[ , !(names(MEs) %in% "MEgrey")]
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.15
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
#need to remove MEgrey colours from dynamis colours
dynamicColors_red <- dynamicColors[dynamicColors!= 'grey']

merge = mergeCloseModules(WGCNA_in3, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
#Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

####
nGenes = ncol(WGCNA_in3)
nSamples = nrow(WGCNA_in3)


# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNA_in3, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, WGCNA_meta_2$wgcna, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 7)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),"(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = "Betty_Module_trait_corr_test.pdf", width = 4, height = 8);
par(mar = c(6, 9, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "I vs NI",
               yLabels = names(MEs),
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.x = 2,
               cex.lab.y = 1,
               cex.text = 1,
               zlim = c(-1,1))
dev.off()

##getting genes associated with module
# names (colors) of the modules
#retrieving genes for MEnavajowhite1, bc of significant correlation 0.75, p=0.05
#20 genes, 15 annotated (Biomart)
#correlates with NI
modNames = substring(names(MEs), 3)
module = "navajowhite1"
column = match(module, modNames)
moduleGenes = moduleColors==module
#Type II diabetes mellitus WP1584 (Wikkipathways human 2019,0.016)
#potassium ion binding (GO:0030955), 0.045

MEnavajowhite1_genes <- data.frame(names(WGCNA_in3)[moduleColors=="navajowhite1"])
write.table(MEnavajowhite1_genes,"MEnavajowhite1_genes.txt")
ME_turquoise_bimart <- read.csv("METurquoise_biomaRt.txt", header = T, stringsAsFactors = F)
METurquoise_genes_annot <- merge(METurquoise_genes,ME_turquoise_bimart,by.x="names.WGCNA_in3..moduleColors.....turquoise..",by.y="Gene.stable.ID.version")
write.table(METurquoise_genes_annot,"METurquoise_genes_annot.txt")

##and retrieving MEdarkseagreen4, -0.77, p-0.04
#correlates with I
module = "darkseagreen4"
column = match(module, modNames)
moduleGenes = moduleColors==module
MEdarkseagreen4_genes <- data.frame(names(WGCNA_in3)[moduleColors=="darkseagreen4"])
write.table(MEdarkseagreen4_genes,"MEdarkseagreen4_genes.txt")

##and retrieving MEindianred4, -0.85, p=0.02
#correlates with I
module = "indianred4"
column = match(module, modNames)
moduleGenes = moduleColors==module
MEindianred4_genes <- data.frame(names(WGCNA_in3)[moduleColors=="indianred4"])
write.table(MEindianred4_genes,"MEindianred4_genes.txt")
