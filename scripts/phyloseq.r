# phyloseq.r
# make global phyloseq object for all miseq runs

# ----- LIBRARIES -----

library(phyloseq)
library(dada2)

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
                "20200228_0094",
                "20201104_0112",
                "20201112_0114")

i = 1
for(run in miseq.runs.emp){
  print(run)
  
  seqtab.filename <- paste("~/Projects/mime-16s/results/", run, "/dada2-merge-chimera-taxo/seqtab_final.rds", sep = "")
  taxa.filename <- paste("~/Projects/mime-16s/results/", run, "/dada2-merge-chimera-taxo/tax_final.rds", sep = "")

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

# Issues with PCR replicates?
sample_data(ps)[is.na(sample_data(ps)$cruise) == TRUE,]

# make a categorial variable for north and south
sample_data(ps)$region <- NA
sample_data(ps)[sample_data(ps)$transect %in% c("LB", "FX", "SB", "IH", "SK"),]$region <- "South"
sample_data(ps)[sample_data(ps)$transect %in% c("KG", "HB", "SI", "MS", "LN", "LA", "KR"),]$region <- "North"


saveRDS(ps, file = "/Users/Clara/Projects/mime-16s/global-ps.rds")

