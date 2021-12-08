#!/usr/bin/env Rscript
# phyloseq.r
# make global phyloseq object for all miseq runs

# ----- LIBRARIES -----

library(phyloseq)
library(dada2)
library(dplyr) # for left_join

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

# run 20190508_0074

miseq.runs.emp <- c("20190503_0073",
                "20190508_0074",
                "20190915_0082",
                "20191002_0084",
                "20191016_0085",
                "20200228_0094",
                "20200331_0098",
                "20200407_0100",
                "20200421_0102",
                "20201104_0112",
                "20201112_0114")

i = 1
for(run in miseq.runs.emp){
  print(run)

  seqtab.filename <- paste("~/Projects/mime-16s/results/", run, "/seqtab_final.rds", sep = "")
  taxa.filename <- paste("~/Projects/mime-16s/results/", run, "/tax_final.rds", sep = "")

  seqtab <- readRDS(seqtab.filename)
  tax_info <- readRDS(taxa.filename)

  # extract info fron sample names
  samples.out <- rownames(seqtab)
  stn <- sapply(strsplit(samples.out, "-"), `[`, 1)
  smp.num <- paste0(sapply(strsplit(samples.out, "-"), `[`, 1), "-", sapply(strsplit(samples.out, "-"), `[`, 2))
  smp.num <- substr(smp.num, 1,9)
  primer <- sapply(strsplit(samples.out, "-"), `[`, 3)
  samdf <- data.frame(stn=stn, smp.num=smp.num, primer=primer, run=run)
  rownames(samdf) <- samples.out

  # build sample info dataframe
  samdf <- left_join(samdf, metadata, copy = FALSE, stringsAsFactors = FALSE)
  rownames(samdf) <- samples.out

  # Issues with PCR replicates?
  samdf[is.na(samdf$cruise) == TRUE,]

  print(i)
  if (i == 1){
    print("first")
    ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(tax_info))
    print(ps)
  } else {
    ps2 <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(tax_info))
    ps <- merge_phyloseq(ps, ps2)
    print(ps)
  }
  i = i + 1
}

# Create new variable for each transect
sample_data(ps)$transect <- substr(sample_data(ps)$stn, 1, 2)

# set dates
sample_data(ps)$date <- as.Date(paste(sample_data(ps)$year,
                                      sample_data(ps)$month,
                                      sample_data(ps)$day, sep = "-"))

# ----- LABEL CONTROLS ----

# positive control = mock
sample_data(ps)[sample_data(ps)$transect == "Mo",]$transect <- "MOCK"
sample_data(ps)[sample_data(ps)$transect == "MOCK",]$stn.num <- "MK000"
sample_data(ps)[sample_data(ps)$transect == "MOCK",]$stn.name <- "Mock community"

#  controls
sample_data(ps)[sample_data(ps)$transect == "Co",]$transect <- "CTRL"
sample_data(ps)[sample_data(ps)$transect == "CTRL",]$stn.num <- "CTRL"
sample_data(ps)[sample_data(ps)$transect == "CTRL",]$stn.name <- "Control"
sample_data(ps)[sample_data(ps)$transect == "CTRL",]$stn.num <- "CTRL"
sample_data(ps)[sample_data(ps)$transect == "CTRL",]$smp.num <-
  substr(rownames(sample_data(ps)[sample_data(ps)$transect == "CTRL",]), 8, 14)


sample_data(ps)$transect <- as.factor(sample_data(ps)$transect)

# Issues with PCR replicates?
sample_data(ps)[is.na(sample_data(ps)$cruise) == TRUE,]

# Information about the primer
sample_data(ps)$primer <- "EMP"

# make a categorial variable for north and south
sample_data(ps)$region <- NA
sample_data(ps)[sample_data(ps)$transect %in% c("LB", "FX", "SB", "IH", "SK", "RS", "HD"),]$region <- "South"
sample_data(ps)[sample_data(ps)$transect %in% c("KG", "HB", "SI", "MS", "LN", "LA", "KR", "EY", "MI", "HF"),]$region <- "North"
sample_data(ps)[is.na(sample_data(ps)$region) == TRUE,]

# variable as factor
sample_data(ps)$cruise <- factor(sample_data(ps)$cruise, levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B7-2017", "B11-2017", "B3-2018", "B7-2018", "B10-2018"))
sample_data(ps)$run <- factor(sample_data(ps)$run)

sample_data(ps)$iscar.nb <- substr(rownames(sample_data(ps)), 1, 9)

# pelagic zone
sample_data(ps)$zone <- cut(sample_data(ps)$depth, breaks = c(-Inf, 0, 200, 1000, Inf), labels = c("surface", "photic", "aphotic", "bathypelgic"))

# ----- SAVE EMP PS -----

saveRDS(ps, file = "/Users/Clara/Projects/mime-16s/global-ps-emp-gtdb.rds")

# ----- ARCHAEA -----

miseq.runs.ark <- c("20200306_0095",
                    "20200318_0097",
                    "20200421_0102",
                    "20200429_0103",
                    "20200506_0104"
                    )

i = 1
for(run in miseq.runs.ark){
  print(run)

  seqtab.filename <- paste("~/Projects/mime-16s/results/", run, "/seqtab_final.rds", sep = "")
  taxa.filename <- paste("~/Projects/mime-16s/results/", run, "/tax_final.rds", sep = "")

  seqtab <- readRDS(seqtab.filename)
  tax_info <- readRDS(taxa.filename)

  # extract info fron sample names
  samples.out <- rownames(seqtab)
  stn <- sapply(strsplit(samples.out, "-"), `[`, 1)
  smp.num <- paste0(sapply(strsplit(samples.out, "-"), `[`, 1), "-", sapply(strsplit(samples.out, "-"), `[`, 2))
  smp.num <- substr(smp.num, 1,9)
  primer <- sapply(strsplit(samples.out, "-"), `[`, 3)
  samdf <- data.frame(stn=stn, smp.num=smp.num, primer=primer, run=run)
  rownames(samdf) <- samples.out

  # build sample info dataframe
  samdf <- left_join(samdf, metadata, copy = FALSE, stringsAsFactors = FALSE)
  rownames(samdf) <- samples.out

  # Issues with PCR replicates?
  samdf[is.na(samdf$cruise) == TRUE,]

  print(i)
  if (i == 1){
    print("first")
    ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(tax_info))
    print(ps)
  } else {
    ps2 <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                    sample_data(samdf),
                    tax_table(tax_info))
    ps <- merge_phyloseq(ps, ps2)
    print(ps)
  }
  i = i + 1
}


sample_data(ps)$transect <- substr(sample_data(ps)$stn, 1, 2)
sample_data(ps)[sample_data(ps)$transect == "Mo",]$transect <- "MOCK"
sample_data(ps)$transect <- as.factor(sample_data(ps)$transect)
sample_data(ps)[sample_data(ps)$transect == "MOCK",]$stn.num <- "MK000"
sample_data(ps)[sample_data(ps)$transect == "MOCK",]$stn.name <- "Mock community"
sample_data(ps)$date <- as.Date(paste(sample_data(ps)$year, sample_data(ps)$month, sample_data(ps)$day, sep = "-"))

# Issues with PCR replicatess?
sample_data(ps)[is.na(sample_data(ps)$cruise) == TRUE,]

sample_data(ps)$primer <- "ARK"

# make a categorial variable for north and south
sample_data(ps)$region <- NA
sample_data(ps)[sample_data(ps)$transect %in% c("LB", "FX", "SB", "IH", "SK", "RS"),]$region <- "South"
sample_data(ps)[sample_data(ps)$transect %in% c("KG", "HB", "SI", "MS", "LN", "LA", "KR"),]$region <- "North"

# variable as factor
sample_data(ps)$cruise <- factor(sample_data(ps)$cruise, levels =c("B8-2010", "B4-2011", "B5-2012", "B3-2013", "B4-2014", "B4-2015","B9-2016", "B11-2017", "B3-2018", "B7-2018"))
sample_data(ps)$run <- factor(sample_data(ps)$run)

sample_data(ps)$iscar.nb <- substr(rownames(sample_data(ps)), 1, 9)

# pelagic zone
sample_data(ps)$zone <- cut(sample_data(ps)$depth, breaks = c(-Inf, 0, 200, 1000, Inf), labels = c("surface", "photic", "aphotic", "bathypelgic"))


# ------- CLEANUP DATASET ----

# ps <- subset_samples(ps, transect != "Re" |
#                       transect != "Ne" |
#                       transect != "Po" |
#                       transect != "Co" |
#                       transect != "Ba")

saveRDS(ps, file = "/Users/Clara/Projects/mime-16s/global-ps-ark-gtdb.rds")
