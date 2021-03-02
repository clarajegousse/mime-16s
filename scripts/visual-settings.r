# visual-settings.r
# To get consistently pretty-looking plots :)

# ----- LIBRARIES -----

library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# ----- COLOURS -----

grapefruit <- "#ED5565"
Grapefruit <- "#DA4453"
bittersweet <- "#FC6E51"
Bittersweet <- "#E9573F" #
orange <- "#EBC9A6"
Orange <- "#DF8830"
sunflower <- "#FFCE54"
Sunflower <- "#F6BB42"
grass <- "#A0D468"
Grass <- "#8CC152"
mint <- "#48CFAD"
Mint <- "#37BC9B"
aqua <- "#4FC1E9"
Aqua <- "#3BAFDA"
jeans <- "#5D9CEC"
Jeans <- "#4A89DC"
lavender <- "#AC92EC"
Lavender <- "#967ADC"
rose <- "#EC87C0"
Rose <- "#D770AD"
lightgrey <- "#F5F7FA"
LightGrey <- "#E6E9ED"
mediumgrey <- "#CCD1D9"
MediumGrey <- "#AAB2BD"
darkgrey <- "#656D78"
DarkGrey <- "#434A54"
chocolate <- "#8D6E63"
Chocolate <- "#6D4C41"

Hydrogen <- "#E6E9ED" # LightGrey
Carbon <- "#434A54" # DarkGrey
Nitrogen <- "#4A89DC" # Jeans
Oxygen <- "#DA4453" # Grapefruit
Fluorine <- "#8CC152" # Grass
Chlorine <- "#8CC152" # Grass
Iodine <- "#967ADC" # Lavender
Phosphorus <- "#FF7F50" # Coral
Sulfur <- "#F6BB42" # Surfur
Iron <- "#E9573F" # Bittersweet
Titanium <- "#AAB2BD" # MediumGrey

# ----- PALETTES -----

Palette1 <- c(Grapefruit, Bittersweet, Orange,
  Sunflower, Grass, Mint, Aqua, Jeans, Lavender, Rose,
  LightGrey, MediumGrey, DarkGrey, Chocolate)

ark.phyla <- c("Aenigmarchaeota", "Asgardarchaeota", "Crenarchaeota", "Euryarchaeota",
                 "Halobacterota", "Hydrothermarchaeota", "Iainarchaeota", "Micrarchaeota",
                 "Nanoarchaeota", "Thermoplasmatota")
ark.phyla.palette <- c(DarkGrey, Grass, Aqua, Mint, Rose, Jeans, Sunflower, Orange, LightGrey, Bittersweet)
names( ark.phyla.palette) <- ark.phyla

# ----- FUNCTIONS -----

# to define colours scales
tax_color_scale <- function(ps, rank) {
  if(rank == "Kingdom"){
    k.colours <- c(Rose, Lavender, Mint, LightGrey)
    names(k.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))
    k.color.scale <- scale_color_manual(name = "Kingdom", values = k.colours)
    return(k.color.scale)
  } else if (length(get_taxa_unique(ps)) == 1 & get_taxa_unique(ps) == "Archaea" & rank == "Phylum") {
    p.colours <- c(DarkGrey, Jeans, Aqua, Grass, MediumGrey, Grapefruit, Sunflower, Mint, Orange)
    names(p.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = rank))
    p.colours.scale <- scale_color_manual(name = rank, values = p.colours)
    return(p.colours.scale)
  } else {
    nb.cols <- 1 + length(get_taxa_unique(ps, taxonomic.rank = rank))
    tax.colours <- colorRampPalette(Palette1[-c(11,12)], bias = 2)(nb.cols)
    names(tax.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = rank))
    tax.color.scale <- scale_color_manual(name = rank, values = tax.colours)
    return(tax.color.scale)
  }
}

# Archaeal Phyla
# Aenigmatarchaeota <- DarkGrey
# Altarchaeota
# Asgardarchaeota <- Jeans
# Hadarchaeota (6)
# Halobacteriota (992)
# Huberarchaeota
# Hydrothermarchaeota <- Aqua
# Iainarchaeota (20)
# Methanobacteriota (337)
# Micrarchaeota (43)
# Nanoarchaeota (207)
# Nanohaloarchaeota (12)
# Thermoplasmatota (770)
# Thermoproteota (579)
# Undinarchaeota (5)
# EX4484-52 (3)
# PWEA01 (6)
# QMZS01 (2)

  tax_fill_scale <- function(ps, rank) {
    if(rank == "Kingdom"){
      k.colours <- c(Rose, Lavender, Mint, LightGrey)
      names(k.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))
      k.fill.scale <- scale_fill_manual(name = "Kingdom", values = k.colours)
      return(k.fill.scale)
    } else {
      nb.cols <- 1 + length(get_taxa_unique(ps, taxonomic.rank = rank))
      tax.colours <- colorRampPalette(Palette1, bias = 1)(nb.cols)
      names(tax.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = rank))
      tax.fill.scale <- scale_fill_manual(name = rank, values = tax.colours)
      return(tax.fill.scale)
    }
  }

# ----- GGPLOT THEME -----

clean_theme = theme_pubr() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        #text = element_text(family = "Muli"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length = unit(0.2,"cm")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 9) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10)


plot_bar2 <-  function (physeq, x, y, fill = NULL, title = NULL, facet_grid = NULL)
{
  mdf = psmelt(physeq)
  p <- ggplot(mdf, aes_string(x = x, y = y, fill = fill)) +
  geom_bar(stat = "identity", position = "stack") + clean_theme +
  tax_fill_scale(physeq, fill) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_line(size = .1, colour = DarkGrey, linetype = 3))
  return(p)
}
