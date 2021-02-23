# phylumanalysis.r


# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

mytheme = theme_pubr()  + 
  theme( #aspect.ratio=1,
        panel.border = element_rect(colour = DarkGrey, fill = NA, size = .75),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.line = element_line(size = 0, color = "red"),
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(color = MediumGrey, face = "italic", size = 8),
        axis.ticks = element_line(size = .75, color = DarkGrey),
        axis.ticks.length = unit(0.2, "cm"),
        strip.background = element_rect(color = "white", size = .75, fill = "white"),
        strip.text.x = element_text(color = DarkGrey, face = "bold"),
        panel.grid.major = element_line(size = .1, colour = MediumGrey),
        text = element_text(family = "Muli")) +
        font("xylab", size = 20, face = "bold") +
        font("xy", size = 12, face = "bold") +
        font("xy.text", size = 9) +
        font("legend.title", size = 12, face = "bold") +
        font("legend.text", size = 10) 

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
library(marmap)

# ----- MAP OF ICELAND -----

xlim <- c(-30, -10) 
ylim <- c(68.5, 62)

# get bathymetry data from NOAA
depth <- getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                       lat1 = ylim[1], lat2 = ylim[2],
                       resolution = 5)

# turn the object into a data.frame
df.depth <- fortify(depth)
iceland <- map_data("world", region = "Iceland")

map <- ggplot() +
  theme_bw() +
  geom_contour(data = depth, aes(x, y, z = z),
               breaks=c(-25, -50, -100, -200, -400),
               colour = MediumGrey, size=0.1) +
  geom_path(data = iceland, aes(long, lat), size=.1, color = DarkGrey) 

# ----- LOAD PHYLOSEQ OBJ -----

ps <- readRDS("/Users/Clara/Projects/mime-16s/global-ps.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

sample_data(ps)$cruise <- factor(sample_data(ps)$cruise, levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B11-2017", "B3-2018", "B7-2018"))
sample_data(ps)$run <- factor(sample_data(ps)$run)
sample_data(ps)$iscar.id <- substr(rownames(sample_data(ps)), 1,9)

df.sample.data <- as.data.frame(sample_data(ps))

# ----- Number of samples per run -----

ggplot(data = df.sample.data , aes(cruise)) +
  geom_bar( fill = DarkGrey) + facet_wrap(~ run, nrow = 11) + 
  mytheme +
  labs(title ="Sample counts per cruise for each MiSeq run",
       subtitle = "with EMP primers",
       caption = "*NA are samples labeled 'Rename' that cannot be associated with sampling metadata.")


nb.cols <- length(levels(df.sample.data$cruise)) + 1
cruise.color <- colorRampPalette(intensePalette)(nb.cols)

map +
  stat_sum(data = df.sample.data[is.na(df.sample.data$cruise) == FALSE,], 
             aes(x = lon, y = lat, size = cruise), alpha = .5) + 
  facet_wrap(~ cruise, ncol = 5) +
  mytheme + theme(legend.position = "bottom") + 
  labs(x = NULL, y = NULL) + guides(size=FALSE) +
  scale_color_manual(values = cruise.color) +
  labs(title ="Number of samples counts per cruise",
       subtitle = "with EMP primers")

map +
  geom_count(data = df.sample.data[is.na(df.sample.data$cruise) == FALSE,], 
           aes(x = lon, y = lat, size = cruise), alpha = .5, color = MediumGrey) + 
  mytheme + theme(legend.position = "bottom") + 
  labs(x = NULL, y = NULL) + guides(size=FALSE) +
  labs(title ="Number of samples counts per cruise",
       subtitle = "with EMP primers")


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

plot_bar(ps3, x="cruise", fill = "Kingdom") +
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
  font("legend.text",size = 10)

