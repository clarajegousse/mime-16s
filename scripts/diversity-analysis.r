# divanalysis.r

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
library(marmap)

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source_url("https://raw.githubusercontent.com/clarajegousse/mime-16s/main/scripts/visual-settings.r")

# ----- MAP OF ICELAND -----

xlim <- c(-30, -10)
ylim <- c(70, 62)

# get bathymetry data from NOAA
depth <- getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                       lat1 = ylim[1], lat2 = ylim[2],
                       resolution = 5)

# turn the object into a data.frame
df.depth <- fortify(depth)
iceland <- map_data("world", region = "Iceland")

# ----- LOAD PHYLOSEQ OBJ -----

ps <- readRDS("/Users/Clara/Projects/mime-16s/global-ps-emp.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# ----- REMOVE THE MOCK -----

ps <- ps %>%
  subset_samples(stn.num != "MK000") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# ----- REMOVE SAMPLES WITH VERY FEW READS -----

good.samples <- rownames(as.data.frame(sample_sums(ps)[sample_sums(ps) > 500 ]))

ps <- ps %>%
  subset_samples( rownames(sample_data(ps)) %in% good.samples) 

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

nb.cols <- 5 + length(get_taxa_unique(ps2, taxonomic.rank = "Phylum"))
phylum.color <- colorRampPalette(Palette1)(nb.cols)

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

ps3 = tax_glom(ps2, taxrank = taxGlomRank, NArm = TRUE)

# ------- SUBSET SPECIFIC SAMPLES ----

# here let's select only the surface samples

ps4 <- subset_samples(ps2, depth =="0")
ps4 <- ps3

# ----- SUBSET SPECIFIC GROUP -----

ps5 <- subset_taxa(ps4, Class %in% c("Gammaproteobacteria"))
ps5

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps5))
standf = function(x, t=total) round(t * (x / sum(x, na.rm = TRUE)))
ps6 = transform_sample_counts(ps5, standf)

plot_bar(ps6, fill = "Order")

sum(is.na(otu_table(ps6)))
psmelt(ps5) %>%
  filter(is.na(Abundance))


# Basic bar graph based on Order

# set the colours
nb.cols <- length(unique(tax_table(ps6)[,"Order"])) + 1
mycolors <- colorRampPalette(Palette1)(nb.cols)

plot_bar(ps6, x="Sample", fill = "Order") +
  scale_fill_manual(values=rev(mycolors)) + 
  clean_theme


# ----- HEATMAP -----

# very cluttered heatmap
plot_heatmap(ps6, method = "NMDS", distance = "bray", na.value = "white")

# For example one can only take ASVs that represent at least 20% of reads in at least one sample
# Remember we normalized all the sampples to median number of reads (total).
# We are left with only 6 ASVs which makes the reading much more easy.

ps7 <- filter_taxa(ps6, function(x) sum(x > total*0.20) > 0, TRUE)
ps7
otu_table(ps7)[1:8, 1:5]

plot_heatmap(ps7, method = "NMDS", distance = "bray")

plot_heatmap(ps7, method = "MDS", distance = "(A+B-2*J)/(A+B-J)",
             taxa.label = "Order", taxa.order = "Order",
             trans=NULL, low="white", high=Jeans, na.value="white") +
    clean_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                 legend.position="right")

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# A+B-2*J “quadratic” squared Euclidean
# A+B-2*J “minimum” Manhattan
# (A+B-2*J)/(A+B) “minimum” Bray-Curtis
# (A+B-2*J)/(A+B) “binary” Sørensen
# (A+B-2*J)/(A+B-J) “binary” Jaccard

# ---- Alpha diversity ----

plot_richness(ps5, measures=c("Chao1", "Shannon")) +
  clean_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position="right")

plot_richness(ps.gamma, measures=c("Chao1", "Shannon"), x = "cruise", color="transect") +
  clean_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  legend.position="right")

# ----- Ordination -----

ps.gamma.ord <- ordinate(ps6, "NMDS", "bray")

plot_ordination(ps6, ps.gamma.ord, type="taxa", color="Order", shape= "year",
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
  mytheme + theme(legend.position = "right") #+ facet_wrap(~cruise)

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

# ----- MAP ----

library(marmap)

xlim <- c(-30, -10)
ylim <- c(70, 62)

# get bathymetry data from NOAA
depth <- getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                       lat1 = ylim[1], lat2 = ylim[2],
                       resolution = 1)

# turn the object into a data.frame
df.depth <- fortify(depth)
iceland <- map_data("world", region = "Iceland")

# select specific taxonomic group
ps.thio <- subset_taxa(ps4, Family %in% c("Thioglobaceae"))

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps.thio))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.thio = transform_sample_counts(ps.thio, standf)

df.thio <- as.data.frame(otu_table(ps.thio))
head(df.thio)

df.thio$sample <-rownames(df.thio)
df.thio$stn <- substr(rownames(df.thio), 1,5)

df.thio$lat <- sample_data(ps.thio)[sample_data(ps.thio)$stn == df.thio$stn,]$lat
df.thio$lon <- sample_data(ps.thio)[sample_data(ps.thio)$stn == df.thio$stn,]$lon
df.thio$cruise <- sample_data(ps.thio)[sample_data(ps.thio)$stn == df.thio$stn,]$cruise

mdf.thio <- melt(df.thio, id=c("sample", "stn","lat", "lon", "cruise"))
head(mdf.thio)

mdf.thio$cruise <- factor(mdf.thio$cruise, levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B11-2017", "B3-2018", "B7-2018"))

xlim <- c(-30, -10)
ylim <- c(62, 70)

ggplot() +
  theme_bw() +
  geom_contour(data = depth, aes(x, y, z = z),
               breaks=c(-25, -50, -100, -200, -400),
               colour="black", size=0.1) +
  geom_path(data = iceland, aes(long, lat), size=0.1) +
  geom_point(data = mdf.thio, aes(x = lon,
               y = lat,
               size = value, color = variable),
             alpha = 0.5) +
  facet_grid(variable~cruise) + theme(aspect.ratio=1) +
  # facet_grid(~cruise) + theme(aspect.ratio=1) +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
labs(x = NULL, y = NULL)

as.character(refseq(ps)["ASV2"])



library(ggrastr)
ggplot() +
  theme_bw() +
  geom_contour(data = depth, aes(x, y, z = z),
               breaks=c(-25, -50, -100, -200, -400),
               colour="black", size=0.1) +
  geom_path(data = iceland, aes(long, lat), size=0.1) +
  rasterise(geom_point(data = mdf.thio[mdf.thio$variable == "ASV2" & mdf.thio$cruise == "B8-2010",], aes(x = lon,
                                  y = lat,
                                  size = value, color = variable),
             alpha = 0.5)) +
   theme(aspect.ratio=1) +
  # facet_grid(~cruise) + theme(aspect.ratio=1) +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = NULL, y = NULL)



data(nw.atlantic)
atl <- as.bathy(nw.atlantic)
library(lattice)
wireframe(unclass(atl), shade = TRUE, aspect = c(1/2, 0.1))

data(nw.atlantic)
atl <- as.bathy(nw.atlantic)


depth2 <- getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                       lat1 = ylim[1], lat2 = ylim[2],
                       resolution = 5)

library(lattice)
wireframe(unclass(depth2), shade = TRUE, aspect = c(1/2, 0.1))



