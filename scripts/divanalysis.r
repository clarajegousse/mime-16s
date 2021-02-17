# divanalysis.r

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

transect.color = scale_color_manual(values = c("FX" = Amethyst, 
                                               "IH" = Aqua,
                                               "KG" = Grass,
                                               "KR" = Bittersweet,
                                               "LA" = Pumpkin, 
                                               "LB" = Orange,
                                               "LN" = Sunflower, 
                                               "MS" = Mint, 
                                               "SB" = Grapefruit,
                                               "SI" = Jeans, 
                                               "SK" = Lavender))

# ----- LIBRARIES -----

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
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

# run 20190508_0074

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

# ---- BUILD PHYLOSEQ OBJECT -----

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(tax_info))

ps74 <- ps

# run 20190915_0082

seqtab.filename <- "~/Projects/mime-16s/results/20190915_0082/dada2-merge-chimera-taxo/seqtab_final.rds"
seqtab <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/20190915_0082/dada2-merge-chimera-taxo/tax_final.rds"
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

# ---- BUILD PHYLOSEQ OBJECT -----

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(tax_info))

ps82 <- ps

# run 20191002_0084

seqtab.filename <- "~/Projects/mime-16s/results/20191002_0084/dada2-merge-chimera-taxo/seqtab_final.rds"
seqtab <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/20191002_0084/dada2-merge-chimera-taxo/tax_final.rds"
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

# ---- BUILD PHYLOSEQ OBJECT -----

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(tax_info))

ps84 <- ps

# ----- MERGE PHYLOSEQ OBJ -----

ps <- merge_phyloseq(ps74, ps82, ps84)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

ps <- ps %>%
  subset_samples(stn.num != "MK000") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# ----- SELECT PROKARYOTES -----

ps <- subset_taxa(ps, Kingdom %in% c("Archaea", "Bacteria"))

# ------- SUBSET SPECIFIC SAMPLES ----

# here let's select only the surface samples

ps1 <- subset_samples(ps, depth =="0")
ps1 


# ----- SSUBSET SPECIFIC GROUP -----

# let's look at the gammaproteobacteria

ps.gamma <- subset_taxa(ps1, Class %in% c("Gammaproteobacteria"))
ps.gamma

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps.gamma))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.gamma = transform_sample_counts(ps.gamma, standf)

# Basic bar graph based on Order

# set the colours
nb.cols <- length(unique(tax_table(ps.gamma)[,"Order"])) + 1
mycolors <- colorRampPalette(flatUI)(nb.cols)

plot_bar(ps.gamma, x="Sample", fill = "Order") +
  scale_fill_manual(values=rev(mycolors)) + theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="right",
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        text=element_text(family="Muli")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) 

# ----- HEATMAP -----

# very cluttered heatmap
plot_heatmap(ps.gamma, method = "NMDS", distance = "bray")

# For example one can only take ASVs that represent at least 20% of reads in at least one sample
# Remember we normalized all the sampples to median number of reads (total). 
# We are left with only 6 ASVs which makes the reading much more easy.

ps.gamma.abund <- filter_taxa(ps.gamma, function(x) sum(x > total*0.20) > 0, TRUE)
ps.gamma.abund
otu_table(ps.gamma.abund)[1:8, 1:5]

plot_heatmap(ps.gamma.abund, method = "NMDS", distance = "bray")

plot_heatmap(ps.gamma.abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Order", taxa.order = "Order", 
             trans=NULL, low="white", high=Jeans, na.value="white") +
  theme_pubr() + mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                 legend.position="right")

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# A+B-2*J “quadratic” squared Euclidean
# A+B-2*J “minimum” Manhattan
# (A+B-2*J)/(A+B) “minimum” Bray-Curtis
# (A+B-2*J)/(A+B) “binary” Sørensen
# (A+B-2*J)/(A+B-J) “binary” Jaccard

# ---- Alpha diversity ----

plot_richness(ps.gamma, measures=c("Chao1", "Shannon")) +
  mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position="right")

plot_richness(ps.gamma, measures=c("Chao1", "Shannon"), x = "cruise", color="transect") +
  mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  legend.position="right") + transect.color


# ----- Ordination -----

ps.gamma.ord <- ordinate(ps.gamma, "NMDS", "bray")

#ps.gamma.ord <- ordinate(ps.gamma, "PCoA", "bray")
#plot_scree(ps.gamma.ord, "Scree plot for Bray/PCoA")

plot_ordination(ps.gamma, ps.gamma.ord, type="taxa", color="Order", shape= "year", 
                title="ASVs")

plot_ordination(ps.gamma, ps.gamma.ord, "samples", color="transect") + 
  geom_point(size=5) + geom_path() + # scale_colour_hue(guide = FALSE) + 
  guides(color=guide_legend(ncol=3)) +
  transect.color +
  mytheme + theme(legend.position = c(0.8, 0.2))
  
# plot_ordination(ps.gamma, ps.gamma.ord, "samples", axes=c(1, 3),
#                color="transect") + geom_line() + geom_point(size=5)

plot_ordination(ps.gamma, ps.gamma.ord, type="taxa", color="Order", 
                title="ASVs", label="Class") + 
  facet_wrap(~Order, 3) + scale_colour_hue(guide = FALSE)


plot_ordination(ps.gamma, ps.gamma.ord, type="samples", color="transect", 
                shape="cruise", title="Samples") + geom_point(size=3)

plot_ordination(ps.gamma, ps.gamma.ord, type="split", color="Order", 
                shape="cruise", title="biplot", label = "station") +  
  geom_point(size=3)


plot_net(ps.gamma, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.9, color="Order", point_label="Genus")

plot_net(ps.gamma.abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Order", point_label="Genus") 
