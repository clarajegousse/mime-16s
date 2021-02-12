#!/usr/bin/env Rscript

# "./scripts/extract_dada2_results.R {input.seqtab} {input.taxo} {output.asv_seq} {output.asv_counts} {output.asv_tax}"

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("At least four arguments must be supplied (one input and three output files).n", call.=FALSE)
}

library(dada2)

# args <- c("results/dada2/seqtab/20190508_0074/seqtab.nochimeras.RDS",
  "results/dada2/taxa/20190508_0074/taxa.RDS",
  "results/dada2/final/20190508_0074/ASVs.fa",
  "results/dada2/final/20190508_0074/ASVs_counts.tsv")

seqtab <- readRDS(args[1])
tax_info <- readRDS(args[2])

# making and writing out a fasta of our final ASV seqs:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, args[3])

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, args[4], sep="\t", quote=F, col.names=NA)
