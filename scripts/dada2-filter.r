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

fns <- list.files(input.path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(input.path, fns), file.path(output.path, fns),
              truncLen = 240, maxEE = 2, truncQ = 10, rm.phix = TRUE,
              compress = TRUE, verbose = TRUE, multithread = TRUE)
