
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
library(devtools)

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
ps4 <- subset_samples(ps3, depth == "0")

# select cruises for which I have MG data
ps4 <- subset_samples(ps4, year != "2010")
ps4 <- subset_samples(ps4, is.na(cruise) != TRUE)
#ps4 <- ps3

# ----- SUBSET SPECIFIC GROUP -----

# select specific taxonomic group
ps.halieaceae <- subset_taxa(ps4, Family %in% c("Halieaceae"))

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps.halieaceae))
standf = function(x, t=total) round(t * (x / sum(x)))
ps.halieaceae= transform_sample_counts(ps.halieaceae, standf)

# set the colours
nb.cols <- length(unique(tax_table(ps.halieaceae)[,"Genus"])) + 1
mycolors <- colorRampPalette(Palette1)(nb.cols)

plot_bar(ps.halieaceae, x="Sample", fill = "Genus") +
  scale_fill_manual(values=rev(mycolors)) + 
  clean_theme + facet_grid(~ region)

# ----- MAP ----

# transform phyloseq to dataframe
df.halieaceae <- as.data.frame(otu_table(ps.halieaceae))
head(df.halieaceae)

df.halieaceae$sample <-rownames(df.halieaceae)
df.halieaceae$stn <- substr(rownames(df.halieaceae), 1,5)

df.halieaceae$lat <- sample_data(ps.halieaceae)[sample_data(ps.halieaceae)$stn == df.halieaceae$stn,]$lat
df.halieaceae$lon <- sample_data(ps.halieaceae)[sample_data(ps.halieaceae)$stn == df.halieaceae$stn,]$lon
df.halieaceae$cruise <- sample_data(ps.halieaceae)[sample_data(ps.halieaceae)$stn == df.halieaceae$stn,]$cruise

mdf.halieaceae <- melt(df.halieaceae, id=c("sample", "stn","lat", "lon", "cruise"),
                       variable.name = "ASV",
                       value.name = "Abundance")
head(mdf.halieaceae)

mdf.halieaceae$cruise <- factor(mdf.halieaceae$cruise, 
                                levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B11-2017", "B3-2018", "B7-2018"))

mdf.halieaceae$Genus <- NA

for (asv in row.names(tax_table(ps.halieaceae))) {
  print(asv)
  mdf.halieaceae[mdf.halieaceae$ASV == asv,]$Genus <- as.character(tax_table(ps.halieaceae)[asv, "Genus"])  
}

xlim <- c(-30, -10)
ylim <- c(62, 70)

ggplot() +
  theme_bw() +
  geom_contour(data = depth, aes(x, y, z = z),
               breaks=c(-25, -50, -100, -200, -400),
               colour="black", size=0.1) +
  geom_path(data = iceland, aes(long, lat), size=0.1) +
  geom_point(data = mdf.halieaceae, aes(x = lon,
                                  y = lat,
                                  size = Abundance, color = Genus),
             alpha = 0.5) +
  facet_grid(Genus ~ cruise) + theme(aspect.ratio=1) +
  # facet_grid(~cruise) + theme(aspect.ratio=1) +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = NULL, y = NULL)

plot.file.name <- paste0(img.path, today, "_distribution_halieaceae.pdf")
ggsave(
  plot.file.name,
  plot = ggplot2::last_plot(),
  path = NULL,
  scale = 1,
  width = 497.42,
  height = 200,
  units = "mm",
  limitsize = FALSE
)

# to visualise sequences
as.character(refseq(ps)["ASV85"])
