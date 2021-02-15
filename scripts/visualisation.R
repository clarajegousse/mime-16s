# visualisation.R

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source("/Users/Clara/Projects/colors/colors.R")
source("/Users/Clara/Projects/colors/colors2.R")

# ----- LIBRARIES -----

library(phyloseq)

# ----- LOAD DATA -----

seqtab.filename <- "~/Projects/mime-16s/results/dada2/seqtab/20190508_0074/seqtab.collapsed.RDS"
seqtab <- readRDS(seqtab.filename)

taxa.filename <- "~/Projects/mime-16s/results/dada2/taxa/20190508_0074/taxa.RDS"
taxa <- readRDS(taxa.filename)
