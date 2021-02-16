# divanalysis.r

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# ----- LIBRARIES -----

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
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
