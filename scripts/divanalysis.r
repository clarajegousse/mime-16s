# divanalysis.r

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# ----- LIBRARIES -----

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(dada2)
library(reshape2)
library(phyloseq)

# ----- LOAD METADATA -----

metadata <- read.csv("/Users/Clara/Projects/mime/data/all-metadata.csv",
                     sep = ";", dec = ".",
                     na.strings = c("NA", ""), strip.white = TRUE,
                     encoding = "utf-8")

# watch out for duplicates!
metadata$smp.num[duplicated(metadata$smp.num)]
metadata <- metadata[!duplicated(metadata$smp.num), ]
rownames(metadata) <- metadata$smp.num

# ----- LOAD DATA -----

seqtab.filename <- "~/Projects/mime-16s/results/20190508_0074/dada2-merge-chimera-taxo/seqtab_final.rds"
seqtab <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/20190508_0074/dada2-merge-chimera-taxo/tax_final.rds"
tax_info <- readRDS(taxa.filename)

# extract info fron sample names
samples.out <- rownames(seqtab)
stn <- sapply(strsplit(samples.out, "-"), `[`, 1)
smp.num <- paste0(sapply(strsplit(samples.out, "-"), `[`, 1), "-", sapply(strsplit(samples.out, "-"), `[`, 2))
smp.num <- substr(smp.num, 1,9)
primer <- sapply(strsplit(samples.out, "-"), `[`, 3)
samdf <- data.frame(stn=stn, smp.num=smp.num, primer=primer)
rownames(samdf) <- samples.out

# build sample info dataframe
samdf <- left_join(samdf, metadata, copy = FALSE, stringsAsFactors = FALSE)
samdf$transect <- substr(samdf$stn, 1, 2)
samdf[samdf$transect == "Mo",]$transect <- "MOCK"
samdf[samdf$transect == "MOCK",]$stn.num <- "MK000"
levels(samdf$stn.name) <- c(levels(samdf$stn.name),"Mock community")
samdf[samdf$stn.num == "MK000",]$stn.name <- "Mock community"
samdf$date <- as.Date(paste(samdf$year, samdf$month, samdf$day, sep = "-"))
rownames(samdf) <- samples.out

# Issues with PCR replicates?
samdf[is.na(samdf$cruise) == TRUE,]

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(tax_info))


ps

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

ps <- ps %>%
  subset_samples(stn.num != "MK000") %>%
  prune_taxa(taxa_sums(.) > 0, .)
