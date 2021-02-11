#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (input and output file).n", call.=FALSE)
}

library(dada2)

seqtab <- readRDS(args[1])

asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, args[2])
