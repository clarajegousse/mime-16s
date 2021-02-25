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

# ----- FUNCTIONS -----

# to define colours scales

tax_fill_scale <- function(ps, rank) {
  if(rank == "Kingdom"){
    k.colours <- c(Rose, Lavender, Mint, LightGrey)
    names(k.colours) <- levels(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))
    k.fill.scale <- scale_fill_manual(name = "Kingdom", values = k.colours)
    return(k.fill.scale)
  } else {
    nb.cols <- 1 + length(get_taxa_unique(ps0, taxonomic.rank = rank))
    tax.colours <- colorRampPalette(Palette, bias = 2)(nb.cols)
    names(tax.colours) <- levels(get_taxa_unique(ps.ark, taxonomic.rank = rank))
    tax.fill.scale <- scale_fill_manual(name = rank, values = tax.colours)
    return(tax.fill.scale)
  }
}

# ----- GGPLOT THEME -----

clean_theme <- theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        axis.line = element_line(size=0,color="red"),
        axis.ticks = element_line(size=.5,color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        text=element_text(family="Muli")) +
  font("xylab",size = 20, face = "bold") +
  font("xy", size=12, face = "bold") +
  font("xy.text", size = 9) +
  font("legend.title",size = 12, face = "bold") +
  font("legend.text",size = 10)
