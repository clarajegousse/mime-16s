# archaea-comparison.r

# ----- LIBRARIES -----

library(devtools)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dada2)
library(phyloseq)

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source_url("https://raw.githubusercontent.com/clarajegousse/mime-16s/main/scripts/visual-settings.r")

# ----- LOAD PHYLOSEQ OBJ -----

ps.emp <- readRDS("/Users/Clara/Projects/mime-16s/global-ps-emp.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps.emp))
names(dna) <- taxa_names(ps.emp)
ps.emp <- merge_phyloseq(ps.emp, dna)
taxa_names(ps.emp) <- paste0("ASV", seq(ntaxa(ps.emp)))
ps.emp

ps.ark <- readRDS("/Users/Clara/Projects/mime-16s/global-ps-ark.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps.ark))
names(dna) <- taxa_names(ps.ark)
ps.ark <- merge_phyloseq(ps.ark, dna)
taxa_names(ps.ark) <- paste0("ASV", seq(ntaxa(ps.ark)))
ps.ark

# quick check
ps.ark1 <- subset_taxa(ps.ark, Kingdom %in% c("Archaea"))
plot_bar(tax_glom(ps.ark1, taxrank = "Phylum"), 
         x="Sample", fill = "Phylum") +
  scale_fill_manual(values = ark.phyla.colours) + 
  clean_theme

ps <- merge_phyloseq(ps.emp, ps.ark)

# ----- SELECT ARCHAEA ONLY -----

# because these were assigned with Silva
ps0 <- subset_taxa(ps, Kingdom %in% c("Archaea"))

# ------- SUBSET SPECIFIC SAMPLES ----

ps1 <- subset_samples(ps0, cruise == "B5-2012")
ps1 

# ----- BASIC BAR PLOT -----

# Basic bar graph based on Order
plot_bar(ps1, x="transect", fill = "Order") +
  tax_fill_scale(ps1, "Order") + clean_theme +
 facet_grid(primer ~ cruise)

# ---- Alpha diversity ----

# to prune OTUs that are not present in any of the samples 
ps2 <- prune_taxa(taxa_sums(ps1) > 0, ps1)

# set the colours
nb.cols <- length(levels(sample_data(ps2)$transect))
mycolors <- colorRampPalette(Palette1, bias = 2)(nb.cols)

plot_richness(ps2, measures=c("Observed", "Chao1", "Shannon", "Simpson"), x = "primer", color = "transect") +
  scale_color_manual(values=mycolors) +
  clean_theme

# ----- SELECT ARCHAEA ONLY -----

# because these were assigned with Silva
ps.nitro <- subset_taxa(ps2, Order %in% c("Nitrosopumilales"))

# ----- BASIC BAR PLOT -----

# Basic bar graph based on Genus

# set the colours
nb.cols <- length(unique(tax_table(ps.nitro)[,"Order"])) + 1
mycolors <- colorRampPalette(Palette)(nb.cols)

plot_bar(ps.nitro, x="transect", fill = "Family") +
  scale_fill_manual(values=mycolors) +
  theme_pubr() +
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
  font("legend.text",size = 10) + facet_wrap(~primer, ncol = 1)


# ----- SELECT ARCHAEA ONLY -----

# because these were assigned with Silva
ps.crenarchaeota <- subset_taxa(ps1, Phylum %in% c("Crenarchaeota"))

# ----- BASIC BAR PLOT -----

# Basic bar graph based on Genus

# set the colours
nb.cols <- length(unique(tax_table(ps.crenarchaeota)[,"Order"])) + 1
mycolors <- colorRampPalette(Palette)(nb.cols)

plot_bar(ps.crenarchaeota, x="transect", fill = "Order") +
  scale_fill_manual(values=mycolors) +
  theme_pubr() +
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
  font("legend.text",size = 10) + facet_wrap(~primer, ncol = 1)

