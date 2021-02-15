#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

args = commandArgs(trailingOnly=TRUE)

# ----- LIBRARY -----

library(dada2); packageVersion("dada2")

# ----- FILE PATHS -----

# input.path <- "/users/home/cat3/projects/mime-16s/results/20190508_0074/cutadapt"
# output.path <- "/users/home/cat3/projects/mime-16s/results/20190508_0074/dada2-filter"

input.path <- args[1]
output.path <- args[2]

# File parsing
pathF <- input.path
pathR <- input.path

filtpathF <- file.path(output.path)
filtpathR <- file.path(output.path)


fastqFs <- sort(list.files(pathF, pattern="R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
