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

# ----- REMOVE THE MOCK -----

ps <- ps %>%
  subset_samples(stn.num != "MK000") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# ----- SELECT PROKARYOTES ONLY -----

# because these were assigned with Silva
ps0 <- subset_taxa(ps, Kingdom %in% c("Archaea", "Bacteria"))

# ----- PREVALENCE FILTERING -----

# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(ps0),
              MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps0)
ps1

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))
ps2

nb.cols <- length(unique(tax_table(ps2)[,2])) + 1
phylum.color <- colorRampPalette(intensePalette)(nb.cols)

ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum) + theme_pubr() + 
  scale_color_manual(values = phylum.color) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ----- AGGLOMERATE CLOSELY RELATED TAXA -----

# Taxonomic agglomeration
# How many genera are present after filtering?
taxGlomRank = "Genus"
length(get_taxa_unique(ps2, taxonomic.rank = taxGlomRank))

ps3 = tax_glom(ps2, taxrank = taxGlomRank)

# ------- SUBSET SPECIFIC SAMPLES ----

# here let's select only the surface samples

ps4 <- subset_samples(ps2, depth =="0")
ps4 

# ----- SSUBSET SPECIFIC GROUP -----

# let's look at the gammaproteobacteria

ps.gamma <- subset_taxa(ps4, Class %in% c("Gammaproteobacteria"))
ps.gamma

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps.gamma))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.gamma = transform_sample_counts(ps.gamma, standf)

# Basic bar graph based on Order

# set the colours
nb.cols <- length(unique(tax_table(ps.gamma)[,"Order"])) + 1
mycolors <- colorRampPalette(intensePalette)(nb.cols)

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
  font("xy.text", size = 9) +
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
                  legend.position="right") 

# ----- Ordination -----

ps.gamma.ord <- ordinate(ps.gamma, "NMDS", "bray")

plot_ordination(ps.gamma, ps.gamma.ord, type="taxa", color="Order", shape= "year", 
                title="ASVs")

# set the colours
nb.cols <- length(unique(sample_data(ps)$transect)) + 1
transect.color <- colorRampPalette(intensePalette)(nb.cols)

plot_ordination(ps.gamma, ps.gamma.ord, "samples", color="transect") + 
  geom_point(size=5) + geom_path() + # scale_colour_hue(guide = FALSE) + 
  guides(color=guide_legend(ncol=3)) +
  scale_color_manual(values = transect.color) +
  mytheme + theme(legend.position = "right")

plot_ordination(ps.gamma, ps.gamma.ord, "samples", color="region",  label = "transect") + 
  geom_point(size=7, alpha=0.2) + geom_path() + # scale_colour_hue(guide = FALSE) + 
  guides(color=guide_legend(ncol=1)) +
  scale_color_manual(values = c(Jeans, Grapefruit)) +
  mytheme + theme(legend.position = "right") + facet_wrap(~cruise)
  
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

# -----
