#!/usr/bin/env Rscript

# ./scripts/extract_dada2_results.R "results/dada2/seqtab/20190508_0074/seqtab.nochimeras.RDS" \
"results/dada2/final/20190508_0074-ASVs.fa" \
"results/dada2/final/20190508_0074-ASVs_counts.tsv"

args = commandArgs(trailingOnly=TRUE)

library(dada2)

# args <- c("results/dada2/seqtab/20190508_0074/seqtab.nochimeras.RDS",
  "results/dada2/final/20190508_0074-ASVs.fa",
  "results/dada2/final/20190508_0074-ASVs_counts.tsv")

seqtab <- readRDS(args[1])
# tax_info <- readRDS(args[2])

# making and writing out a fasta of our final ASV seqs:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = args[2])

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, args[3], sep="\t", quote=F, col.names=NA)
