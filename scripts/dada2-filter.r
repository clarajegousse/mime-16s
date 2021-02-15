#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

library(dada2); packageVersion("dada2")
# Filename parsing
path <- "/users/home/cat3/projects/mime-16s/results/20190508_0074/cutadapt"
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions
# Filtering
filterAndTrim(file.path(path,fns), file.path(filtpath,fns),
              truncLen=c(240, 200), maxEE=2, truncQ=10, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)