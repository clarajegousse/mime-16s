#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

# args = commandArgs(trailingOnly=TRUE)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Run DADA2 filter")

# Add command line arguments
p <- add_argument(p, "--input_path", help="Input path", type = "character")
p <- add_argument(p, "--output_path", help="Output path", type = "character")
p <- add_argument(p, "--fwd_trunc_len", help="fwd_trunc_len",  type="numeric", default=250)
p <- add_argument(p, "--rev_trunc_len", help="rev_trunc_len",  type="numeric", default=200)
p <- add_argument(p, "--maxee", help="MaxEE",  type="numeric", default=2)
p <- add_argument(p, "--truncq", help="truncQ",  type="numeric", default=11)

# Parse the command line arguments
argv <- parse_args(p)
# print(argv)

# ----- LIBRARY -----

library(dada2); packageVersion("dada2")

# ----- FILE PATHS -----

# File parsing
pathF <- argv$input_path
pathR <- argv$input_path

# pathF <- "/users/home/cat3/projects/mime-16s/results/20200416_0101/cutadapt"
pathR <- "/users/home/cat3/projects/mime-16s/results/20200416_0101/cutadapt"

filtpathF <- file.path(argv$output_path)
filtpathR <- file.path(argv$output_path)

# filtpathF <- "/users/home/cat3/projects/mime-16s/results/20200416_0101/dada2-filter"
# filtpathR <- "/users/home/cat3/projects/mime-16s/results/20200416_0101/dada2-filter"


fastqFs <- sort(list.files(pathF, pattern="R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd = file.path(pathF, fastqFs),
  filt = file.path(filtpathF, fastqFs),
  rev = file.path(pathR, fastqRs),
  filt.rev = file.path(filtpathR, fastqRs),
  truncLen = c(argv$fwd_trunc_len, argv$rev_trunc_len),
  maxEE = argv$maxee, truncQ = argv$truncq,
  compress = TRUE, verbose = TRUE, multithread = TRUE)


# filterAndTrim(fwd = file.path(pathF, fastqFs),
#     filt = file.path(filtpathF, fastqFs),
#     rev = file.path(pathR, fastqRs),
#     filt.rev = file.path(filtpathR, fastqRs),
#     truncLen = c(240,200),
#     #maxEE = 2, truncQ = 2,
#     compress = TRUE, verbose = TRUE, multithread = TRUE)
