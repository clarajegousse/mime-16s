#!/usr/bin/env Rscript

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (input and output file).n", call.=FALSE)
}
args = commandArgs(trailingOnly=TRUE)

library(dada2)

seqtab <- readRDS(args[1])
uniquesToFasta(seqtab, args[2])
