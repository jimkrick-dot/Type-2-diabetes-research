##### SSgSEA Analysis ####
# Load required libraries
library(tidyverse)
library(data.table)
library(GSVA)
# 1.2 Prepare cell markers
cellMarker <- data.table::fread("cellMarker.csv", data.table = F)  # Read cell marker file
colnames(cellMarker)[2] <- "celltype"  # Rename the second column to "celltype"
# Split the marker file into a list grouped by cell type
type <- split(cellMarker, cellMarker$celltype)
# Extract and combine unique gene markers for each cell type
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
# Save the processed cell marker data
save(cellMarker, file = "cellMarker_ssGSEA.Rdata")
## 1.3 Prepare the expression matrix
# Ensure the matrix has genes as rows and samples as columns
expr <- data.table::fread("tpms.txt", data.table = F)  # Read expression file
rownames(expr) <- expr[,1]  # Assign the first column as row names
expr <- expr[,-1]  # Remove the first column
expr <- as.matrix(expr)  # Convert to matrix format
# 2. Perform ssGSEA for immune infiltration quantification
ssgsea_param <- ssgseaParam(
  exprData = expr,   # Input expression matrix
  geneSets = cellMarker,  # Predefined cell marker gene sets
  minSize = 1,  # Minimum number of genes in a gene set
  maxSize = Inf,  # No upper limit on gene set size
  alpha = 0.25,  # Weight parameter for ssGSEA
  normalize = TRUE  # Normalize enrichment scores
)
# Run GSVA with the ssGSEA parameters
gsvadata <- gsva(
  param = ssgsea_param,  # ssGSEA parameter object
  verbose = TRUE  # Print detailed information
)
# Format and save the output
a <- gsvadata %>% t() %>% as.data.frame()
write.table(a, "ssGSEA.txt", sep = "\t", row.names = T, col.names = NA, quote = F)