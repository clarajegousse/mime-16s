#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

# args = commandArgs(trailingOnly=TRUE)
library(argparser, quietly=TRUE)

library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Round a floating point number")

# input_path <- args[1]
# output_path <- args[2]
# trunc.len.fwd <- as.numeric(args[3])
# trunc.len.rev <- as.numeric(args[4])
# maxee <- as.numeric(args[5])
# truncq <- as.numeric(args[6])

# Create a parser
p <- arg_parser("Run DADA2 filter")

# Add command line arguments
p <- add_argument(p, "--input_path", help="Input path", default=FALSE)
p <- add_argument(p, "--output_path", help="Output path", default=FALSE)
p <- add_argument(p, "--trunc_len_fwd", help="number of decimal places",  type="numeric", default=250)
p <- add_argument(p, "--trunc_len_rev", help="number of decimal places",  type="numeric", default=200)
p <- add_argument(p, "--maxee", help="MaxEE",  type="numeric", default=2)
p <- add_argument(p, "--truncq", help="truncQ",  type="numeric", default=11)

# Parse the command line arguments
argv <- parse_args(p)

# ----- LIBRARY -----

library(dada2); packageVersion("dada2")

# ----- FILE PATHS -----

# input.path <- "/users/home/cat3/projects/mime-16s/results/20190508_0074/cutadapt"
# output.path <- "/users/home/cat3/projects/mime-16s/results/20190508_0074/dada2-filter"

# input.path <- args[1]
# output.path <- args[2]
# trunc.len.fwd <- as.numeric(args[3])
# trunc.len.rev <- as.numeric(args[4])
# maxee <- as.numeric(args[5])
# truncq <- as.numeric(args[6])

# File parsing
pathF <- argv$input_path
pathR <- argv$input_path

filtpathF <- file.path(argv$output_path)
filtpathR <- file.path(argv$output_path)

fastqFs <- sort(list.files(pathF, pattern="R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd = file.path(pathF, fastqFs),
  filt = file.path(filtpathF, fastqFs),
  rev = file.path(pathR, fastqRs),
  filt.rev = file.path(filtpathR, fastqRs),
  truncLen = c(argv$trunc_len_fwd, argv$trunc_len_rev),
  maxEE = argv$maxee, truncQ = argv$truncq
  compress = TRUE, verbose = TRUE, multithread = TRUE)
