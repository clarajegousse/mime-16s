# divanalysis.r

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# ----- LIBRARIES -----

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

# ----- LOAD DATA -----

seqtab.filename <- "~/Projects/mime-16s/results/20190508_0074/dada2-merge-chimera-taxo/seqtab_final.rds"
seqtab <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/20190508_0074/dada2-merge-chimera-taxo/tax_final.rds"
tax_info <- readRDS(taxa.filename)

# -----

ps <- phyloseq(tax_table(tax_info),
               otu_table(seqtab, taxa_are_rows = TRUE))


# ----- METADATA -----
# Generate the metadata table from file names

metadata <- read.csv("/Users/Clara/Projects/mime/data/all-metadata.csv",
                     sep = ";", dec = ".",
                     na.strings = c("NA", ""), strip.white = TRUE,
                     encoding = "utf-8")

# watch out for duplicates!
metadata$smp.num[duplicated(metadata$smp.num)]
metadata <- metadata[!duplicated(metadata$smp.num), ]
rownames(metadata) <- metadata$smp.num
