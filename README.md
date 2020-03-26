# Cell Invasion in Digital Microfluidic Microgel Systems
Using digital microfluidics, complex matrices mimicing basement membranes that breast cancer cells invade into were made and breast cancer cells invading these matrices were studied first by using immunohistochemistry and then by performing RNA-seq on specific fractions relating to the breast cancer cells invading into the matrix and the cells residing at the matrix interface.  This is a study of breast cancer cell invasion that probes transcriptomes of invading versus non-invading cells in a 3D matrix created on a digital microfluidic device.

## General Workflow
This RNA-seq pipeline works through normalizing then comparing invading vs non-invading cell populations extracted using this wet lab pipeline (Li et al, 2020) and includes the software pipeline used to look for differentially expressed genes using edgeR and gene correlations using WGCNA.  Other data visualizations include dimensionality reduction using UMAP and T-SNE and hierarchial cluster and visualizations using heatmap.2 from the gplots package.
The following pipelines can be summarized as:
1. Implementing EdgeR: normalizing data to counts per million(CPM), pulling out differentially expressed genes (FDR<0.05), plotting with heatmap.2 package.
2. Implementing WGCNA: variance stabilizing gene matrix, developing thresholds for WGCNA, plotting with ggplot2.
3. Generating statistics (in Supplementary Figures): Gathering and plotting genes versus reads data for each sample, utilizing UMAP package, and summarizing controls versus CIMMS derived RNA-seq data and plotting with heatmap.2 package.

### Publication
Submitted:  
Li, B.B., Scott, E.Y., Chamberlain, M.D., Duong, B.T.V., Zhang, S., Done, S.J. & A.R. Wheeler.  "Cell Invasion in Digital Microfluidic Microgel Systems". 2020

