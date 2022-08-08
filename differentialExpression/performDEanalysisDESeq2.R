# rmf 2.12.2021, last modified 5.11.2022

### ARGUMENTS ###
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  print("Usage: Rscript performDEanalysisDESeq2.R <table of counts> <design table> <outbase>");
  q();
}

countsTable <- args[1]
designTable <- args[2]
outbase <- args[3]

### load libraries ###
library(DESeq2, quietly = T)
library(dplyr, quietly = T)
library(tidyverse, quietly = T)

### read in files ###
countsMatrix <- as.matrix(read.table(countsTable, header = TRUE, row.names = 1, check.names = FALSE))
colnames(countsMatrix)

design <- read.table(designTable, row.names = 1)
head(design)

### wrangle design table ###
colnames(design) <- c("condition")
design <- as.data.frame(design)
head(design)

### round data ###
roundedCounts <- round(countsMatrix)

### create DESeq data set object ###
dds <- DESeqDataSetFromMatrix(countData = roundedCounts, colData = design, design = ~ condition)
print(nrow(dds))

### filter by count ###
# first define boolean vectors which will subset the dataframes into young and old
young <- design[, 1] == "young"
old <- design[, 1] == "old"

# remove miRs with median expression under 1 for either age (under 1 ok for one but not both)
# NOTE: that this filtering is on unnormalized counts
dds_filtered <- dds[rowMedians(counts(dds), cols = young) >= 1 | rowMedians(counts(dds), cols = old) >= 1]  # keep rows with median >= 1
nrow(dds_filtered)

### DEG analysis ###
dds_object <- DESeq(dds_filtered)
sizeFactors(dds_object)

# alpha = 0.05
res0.05 <- results(dds_object, alpha = 0.05, contrast = c("condition", "old", "young"))
summary(res0.05)

res0.05_df <- as.data.frame(res0.05)
mirIDs <- rownames(res0.05_df)
toWrite <- cbind(mirIDs, res0.05_df)

write.table(toWrite, file = paste("DEanalysisResults_alpha0.05", outbase, ".txt", sep=''), sep = '\t', quote = FALSE, row.names = FALSE)

## alpha = 0.1
#res0.1 <- results(dds_object, alpha = 0.1)
#summary(res0.1)

#res0.1_df <- as.data.frame(res0.1)
#mirIDs <- rownames(res0.1_df)
#toWrite <- cbind(mirIDs, res0.1_df)

#write.table(toWrite, file = paste("DEanalysisResults_alpha0.1_", outbase, sep = ''), sep = '\t', quote = FALSE)
