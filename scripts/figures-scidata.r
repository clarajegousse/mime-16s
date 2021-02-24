# figure-scidata.r

# ----- LIBRARIES -----

library(devtools)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(dada2)
library(phyloseq)
library(marmap)

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source_url("https://raw.githubusercontent.com/clarajegousse/colors/master/colors.R")

# ----- MAP OF ICELAND -----

xlim <- c(-30, -5) 
ylim <- c(69, 62)

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
               #breaks=c(-25, -50, -100, -200, -400),
               colour = MediumGrey, size=0.1) +
  geom_path(data = iceland, aes(long, lat), size=.1, color = DarkGrey)
map 

# ----- LOAD PHYLOSEQ OBJ -----

ps.emp <- readRDS("/Users/Clara/Projects/mime-16s/global-ps-emp.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps.emp))
names(dna) <- taxa_names(ps.emp)
ps.emp <- merge_phyloseq(ps.emp, dna)
taxa_names(ps.emp) <- paste0("ASV", seq(ntaxa(ps.emp)))
ps.emp

ps <- ps.emp

# ---- SAMPLING INFO -----

df.sample.data <- as.data.frame(sample_data(ps))

nb.cols <- length(levels(df.sample.data$cruise))
cruise.color <- colorRampPalette(Palette)(nb.cols)

df <- df.sample.data[is.na(df.sample.data$cruise) == FALSE,]

df1 <- unique(df[,c("lat", "lon", 'stn', "cruise")])

fig1a <- map +
  geom_point(data = df,
             aes(x = lon, y = lat, color = cruise),
             shape = 1) + 
  geom_text(data = df, aes(x = lon, y = lat, label = stn), 
            size = 1, color = DarkGrey) +
  facet_wrap(~ cruise, ncol = 5) +
  labs(x = NULL, y = NULL) + guides(size=FALSE) +
  scale_color_manual(values = cruise.color) +
  labs(title ="Sampling stations over the years") +
  theme_pubr()  + 
  theme( #aspect.ratio=1,
    panel.border = element_rect(colour = DarkGrey, fill = NA, size = .75),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line = element_line(size = 0, color = "red"),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(color = MediumGrey, face = "italic", size = 8),
    axis.ticks = element_line(size = .1, color = DarkGrey),
    axis.ticks.length = unit(0.1, "cm"),
    strip.background = element_rect(color = "white", size = .75, fill = "white"),
    strip.text.x = element_text(color = DarkGrey, face = "bold"),
    panel.grid.major = element_line(size = .1, colour = MediumGrey),
    text = element_text(family = "Muli")) +
  font("xylab", size = 20, face = "bold") +
  font("xy", size = 12, face = "bold") +
  font("xy.text", size = 5) +
  font("legend.title", size = 12, face = "bold") +
  font("legend.text", size = 10) 

fig1a
# 

df2 <- df[,c("lat", "lon", 'stn', 'transect', 'cruise', 'depth')]
df2$stn.num <- substr(df2$stn, 5, 5)

ggplot(data = df2) +
  geom_point(aes(x = stn.num, y = depth, color = cruise)) +
  facet_grid(transect ~ cruise) + 
  scale_color_manual(values = cruise.color) +
  scale_y_continuous(trans = "reverse") +
  theme_pubr()  + 
  theme( #aspect.ratio=1,
    panel.border = element_rect(colour = DarkGrey, fill = NA, size = .75),
    axis.text.x = element_text(vjust = 0.5, hjust = 1),
    axis.line = element_line(size = 0, color = "red"),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(color = MediumGrey, face = "italic", size = 8),
    axis.ticks = element_line(size = .1, color = DarkGrey),
    axis.ticks.length = unit(0.1, "cm"),
    strip.background = element_rect(color = "white", size = .75, fill = "white"),
    strip.text.x = element_text(color = DarkGrey, face = "bold"),
    panel.grid.major = element_line(size = .1, colour = MediumGrey),
    text = element_text(family = "Muli")) +
  font("xylab", size = 20, face = "bold") +
  font("xy", size = 12, face = "bold") +
  font("xy.text", size = 5) +
  font("legend.title", size = 12, face = "bold") +
  font("legend.text", size = 10) 

