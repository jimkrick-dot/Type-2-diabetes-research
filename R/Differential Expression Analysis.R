##### Differential Expression Analysis (DEG) #####
library(DESeq2)
library(tidyverse)
# Read gene count matrix
counts <- read.table("counts.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
# Define group information (first 50 columns as "T", others as "N")
conditions <- data.frame(
  sample = colnames(counts),
  group = factor(ifelse(seq_len(ncol(counts)) <= 50, "T", "N"), levels = c("N", "T"))
) %>% column_to_rownames("sample")
# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conditions, design = ~ group)
# Run differential expression analysis
dds <- DESeq(dds)
# Extract DEG results
res <- results(dds)  # Set the adjusted p-value threshold
save(res, file = "DEG.Rda")
