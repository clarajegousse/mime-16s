#!/usr/bin/env Rscript
# dada2-filter.R input_path output_path trunc_len maxee truncq

# args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-i", "--input_path", action="store_true", default=FALSE,
  help="input path [default: none]")
parser$add_argument("-o", "--output_path", action="store_true", default=FALSE,
  help="output path [default: none]")
parser$add_argument("-f", "--trunclenfwd", type="integer", action="store_true", default=250,
  help="truncLen forward [default %(default)s]")
parser$add_argument("-r", "--trunclenrev", type="integer", action="store_true", default=200,
  help="truncLen reverse [default %(default)s]")
parser$add_argument("-m", "--maxee", type="integer", action="store_true", default=2,
    help="MaxEE [default %(default)s]")
parser$add_argument("-t", "--truncq", type="integer", action="store_true", default=11,
        help="TruncQ [default %(default)s]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

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
pathF <- input_path
pathR <- input_path

filtpathF <- file.path(output_path)
filtpathR <- file.path(output_path)

fastqFs <- sort(list.files(pathF, pattern="R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd = file.path(pathF, fastqFs),
  filt = file.path(filtpathF, fastqFs),
  rev = file.path(pathR, fastqRs),
  filt.rev = file.path(filtpathR, fastqRs),
  truncLen = c(trunc.len.fwd, trunc.len.rev),
  maxEE = maxee, truncQ = truncq
  compress = TRUE, verbose = TRUE, multithread = TRUE)
