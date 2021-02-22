# phylumanalysis.r


# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

mytheme = theme_pubr()  + 
  theme(aspect.ratio=1,
        panel.border = element_rect(colour = "black", fill=NA, size=.8),
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        text=element_text(family="Muli")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) 


# ----- LIBRARIES -----

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(vegan)
library(dplyr)
library(scasles)
library(grid)
library(dada2)
library(reshape2)
library(phyloseq)

# ----- LOAD PHYLOSEQ OBJ -----

ps <- readRDS("/Users/Clara/Projects/mime-16s/global-ps.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

sample_data(ps)$cruise <- factor(sample_data(ps)$cruise, levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B11-2017", "B3-2018", "B7-2018"))

# ----- REMOVE THE MOCK -----

ps <- ps %>%
  subset_samples(stn.num != "MK000") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# ----- SELECT PROKARYOTES ONLY -----

# because these were assigned with Silva
ps0 <- subset_taxa(ps, Kingdom %in% c("Archaea", "Bacteria"))


# ----- AGGLOMERATE CLOSELY RELATED TAXA -----

# Taxonomic agglomeration
# How many genera are present after filtering?
taxGlomRank = "Phylum"
length(get_taxa_unique(ps0, taxonomic.rank = taxGlomRank))

ps3 = tax_glom(ps0, taxrank = taxGlomRank)

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps3))
standf = function(x, t=total) round(t * (x / sum(x)))
ps3 = transform_sample_counts(ps3, standf)

# Basic bar graph based on Order

# set the colours
nb.cols <- length(unique(tax_table(ps3)[,"Phylum"])) + 1
mycolors <- colorRampPalette(intensePalette)(nb.cols)

plot_bar(ps3, x="cruise", fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=rev(mycolors)) +
  scale_color_manual(values=rev(mycolors)) + theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="right",
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        text=element_text(family="Muli")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 9) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) + facet_grid(~Kingdom)

