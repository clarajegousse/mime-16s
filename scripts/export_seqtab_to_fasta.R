#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (input and output file).n", call.=FALSE)
}

library(dada2)

# get sample id from the filename (first argument)
sample <- sapply(strsplit(sapply(strsplit(file, "[/]"), tail, 1), "[.]"), `[`, 1)

seqtab <- readRDS(args[1])

asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste0(">", sample, "-ASV-", i)
}

uniquesToFasta(seqtab, args[2], ids = asv_headers)
