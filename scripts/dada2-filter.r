#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

library(dada2); packageVersion("dada2")
# Filename parsing
path <- "path/to/FWD"
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns),
              truncLen=240, maxEE=1, truncQ=11, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
