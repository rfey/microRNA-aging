# rmf 2.20.2019, last modified 1.14.2022

# creates a string vector with all arguments
args = commandArgs(trailingOnly=TRUE)

usage <- 'arguments: <counts table> <design table> <outbase>'

if (length(args)!= 3){
   stop(usage)
}

countsTable <- args[1]
designTable <- args[2]
outbase <- args[3]

library(DESeq2)
library(edgeR)
library(dplyr, quietly = T)
library(tidyverse, quietly = T)

# read in files
rawCounts <- as.matrix(read.delim(countsTable,row.names=1))
design <- read.table(designTable, row.names = 1)

# wrangle design table
colnames(design) <- c("condition")
design <- as.data.frame(design)

# filtering on counts to be consistent with DESeq method
# first define boolean vectors which will subset the dataframes into young and old
young <- design[, 1] == "young"
old <- design[, 1] == "old"

# remove miRs with median expression under 1 for either age (under 1 ok for one but not both)
roundedCounts <- round(rawCounts)
print(paste('Number of miRs before filtering:',nrow(roundedCounts),sep=" "))
filteredCounts <- roundedCounts[rowMedians(roundedCounts, cols = young) >= 1 | rowMedians(roundedCounts, cols = old) >= 1, ]
print(paste('Number of miRs after filtering:',nrow(filteredCounts),sep=" "))

# create DGE list object
filteredCountsDGE <- DGEList(counts = filteredCounts, group = design$condition)

# calculate normalization factors with TMM norm method
filteredCountsDGE <- calcNormFactors(filteredCountsDGE, method = 'TMM')

# estimate dispersion
counts.disp <- estimateDisp(filteredCountsDGE)

# exact test
et <- exactTest(counts.disp)

# multiple test correction
FDR <- p.adjust(et$table$PValue, method = "BH")
et$table <- cbind(et$table, FDR)

# write outfile
mirID <- row.names(et$table)
toWrite <- cbind(mirID, et$table)

outfile <- paste('resultsDEanalysis',outbase,'.txt', sep='')
write.table(toWrite, file = outfile, sep = '\t', quote = FALSE, row.names = FALSE)
