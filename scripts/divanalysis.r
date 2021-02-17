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
        legend.position = c(0.2, 0.7),
        text=element_text(family="Muli")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 12) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10) 

# ----- LIBRARIES -----

library(ggplot2)
library(ggpubr)
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

# ---- BUILD PHYLOSEQ OBJECT -----

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


par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(taxa_sums(ps), TRUE)[1:N]/nsamples(ps), las=2, names.arg = "Order")

# ------- 

head(tax_table(ps))

ps.gamma <- subset_taxa(ps, Class %in% c("Gammaproteobacteria"))
ps.gamma

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps.gamma))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.gamma = transform_sample_counts(ps.gamma, standf)

plot_bar(ps.gamma, x="Sample", fill = "Order") + facet_wrap(~ year)
plot_bar(ps.gamma, x="stn.name", fill = "Order")


plot_bar(ps.gamma, x="stn.name", fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

plot_heatmap(ps.gamma, method = "NMDS", distance = "bray")

ps.gamma.abund <- filter_taxa(ps.gamma, function(x) sum(x > total*0.20) > 0, TRUE)
ps.gamma.abund
otu_table(ps.gamma.abund)[1:8, 1:5]

plot_heatmap(ps.gamma.abund, method = "NMDS", distance = "bray")

plot_heatmap(ps.gamma.abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Order", taxa.order = "Order", 
             trans=NULL, low=Clouds, high=Grass, na.value="beige")

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# A+B-2*J “quadratic” squared Euclidean
# A+B-2*J “minimum” Manhattan
# (A+B-2*J)/(A+B) “minimum” Bray-Curtis
# (A+B-2*J)/(A+B) “binary” Sørensen
# (A+B-2*J)/(A+B-J) “binary” Jaccard

plot_richness(ps.gamma, measures=c("Chao1", "Shannon"))
plot_richness(ps.gamma, measures=c("Chao1", "Shannon"), x="season", color="transect")

ps.gamma.ord <- ordinate(ps.gamma, "NMDS", "bray")

#ps.gamma.ord <- ordinate(ps.gamma, "PCoA", "bray")
#plot_scree(ps.gamma.ord, "Scree plot for Bray/PCoA")

plot_ordination(ps.gamma, ps.gamma.ord, type="taxa", color="Order", shape= "year", 
                title="ASVs")

plot_ordination(ps.gamma, ps.gamma.ord, "samples", color="transect") + 
  geom_point(size=5) + geom_path() + # scale_colour_hue(guide = FALSE) + 
  mytheme
  

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
         maxdist = 0.7, color="Order", point_label="Genus")

plot_net(ps.gamma.abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Order", point_label="Genus") 
