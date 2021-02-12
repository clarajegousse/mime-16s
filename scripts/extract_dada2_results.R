#!/usr/bin/env Rscript

# ./scripts/export_seqtab_to_fasta.R \
# results/dada2/merged/20190508_0074/FX003-016-16S-V4_S58.RDS
# results/dada2/merged/20190508_0074/FX003-016-16S-V4_S58.fasta

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least two arguments must be supplied (input and output file).n", call.=FALSE)
}

library(dada2)

# args <- c("results/dada2/merged/20190508_0074/FX003-016-16S-V4_S58.RDS", "results/dada2/merged/20190508_0074/FX003-016-16S-V4_S58.fasta")

seqtab <- readRDS(args[1])


# making and writing out a fasta of our final ASV seqs:
asv_seqs <- seqtab$sequence
asv_headers <- vector(dim(seqtab)[1], mode="character")

for (i in 1:dim(seqtab)[1]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, args[2])


  # count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)
