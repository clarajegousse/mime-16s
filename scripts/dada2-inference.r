#!/usr/bin/env Rscript
# dada2-inference.R

# ----- LIBRARIES -----

library(argparser, quietly=TRUE)
library(dada2); packageVersion("dada2")

# ----- ARGUMENT PARSING -----

# Create a parser
p <- arg_parser("Run DADA2 inference")

# Add command line arguments
p <- add_argument(p, "--input_path", help="Input path", type = "character")
p <- add_argument(p, "--output_path", help="Output path", type = "character")

# Parse the command line arguments
argv <- parse_args(p)
# print(argv)

# ----- INFERENCE -----

# File parsing
filtpathF <- argv$input_path
filtpathR <- argv$input_path

filtFs <- list.files(filtpathF, pattern="R1.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R2.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

ifelse(!dir.exists(file.path(argv$output_path)), dir.create(file.path(argv$output_path)), FALSE)
output.filename <- paste(argv$output_path, "/seqtab.rds", sep = "", collapse = NULL)
saveRDS(seqtab, output.filename) # CHANGE ME to where you want sequence table saved
