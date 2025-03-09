library("tidyverse")
library("WGCNA")            
library(DESeq2)

# Load expression matrix
exp <- read.table("tpms_T2D.txt", sep = "\t", row.names = 1, check.names = F, header = T)
load("DEG.rda")

# Process the input expression matrix
DEG <- as.data.frame(res) %>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0.3, pvalue < 0.05)
input <- exp[rownames(DEG), ]
datExpr0 = as.data.frame(t(input))

# Start WGCNA
# Check for missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# If there are issues, filter out problematic genes and samples
if (!gsg$allOK){
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot sample clustering
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 80000, col = "red") # Cut-off line for outlier removal

# Remove samples below the threshold
clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
table(clust)
keepSamples = (clust == 1)
datExpr0 = datExpr0[keepSamples, ]
dev.off()

# Re-cluster samples after filtering
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)

# Count the number of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
save(datExpr0, nGenes, nSamples, file = "Step01-WGCNA_input.Rda")

# Network construction and module identification
enableWGCNAThreads()  # Enable multi-threading
powers = c(1:20)      # Set power values
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot power selection results
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red") # Adjust as needed

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# Construct adjacency matrix
softPower = sft$powerEstimate
softPower = 12
adjacency = adjacency(datExpr0, power = softPower)

# Convert adjacency matrix to TOM matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
save(TOM, file = "TOM.Rda")

# Gene clustering
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# Identify dynamic modules
minModuleSize = 80
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Assign colors to modules
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot gene dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge similar modules
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.0001  # Cut height for module merging
abline(h = MEDissThres, col = "red")

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
dev.off()

# Load and process clinical data
clinical <- read.table("ssGSEAwgcna.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
clinical <- clinical[rownames(datExpr0), ]
identical(rownames(clinical), rownames(datExpr0))

# Convert clinical data to numeric
datTraits = as.data.frame(do.call(cbind, lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# Cluster samples and plot heatmap
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Module-trait correlation analysis
MEs = orderMEs(MEs)
moduleTraitCor = cor(MEs, datTraits, use = "p")
write.table(file="Step04-modPhysiological.cor.xls", moduleTraitCor, sep="\t", quote=F)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.table(file="Step04-modPhysiological.p.xls", moduleTraitPvalue, sep="\t", quote=F)

# Visualize correlation matrix
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()
