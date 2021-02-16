# visualisation.R

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# ----- LIBRARIES -----

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling

# ----- LOAD DATA -----

seqtab.filename <- "~/Projects/mime-16s/results/dada2/seqtab/20190508_0074/seqtab.collapsed.RDS"
seqtab.nochim <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/dada2/taxa/20190508_0074/taxa.RDS"
tax_info <- readRDS(taxa.filename)

# ----- 

ps <- phyloseq(tax_table(tax_info),
               otu_table(seqtab.nochim, taxa_are_rows = TRUE))


asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
head(asv_fasta)
#write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)
# write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- colnames(tax_info)

colnames(tax_info)

# samples names
head(rownames(seqtab.nochim))

# avs names
head(rownames(asv_tab))

# sequences
head(colnames(seqtab.nochim))

# ranks
head(colnames(tax_info))

# sequences
head(rownames(tax_info))

merge(x, y, by, by.x, by.y, sort = TRUE)

match(colnames(tax_info), rownames(t(seqtab.nochim))
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)




