#!/usr/bin/env Rscript
# scripts/dada2-merge-chimera-taxo.r

# ----- LIBRARIES -----

library(dada2); packageVersion("dada2")
library(argparser, quietly=TRUE)

# ----- ARGUMENT PARSING -----

# Create a parser
p <- arg_parser("Run DADA2 Merge, remove chimera and assign taxo")

# Add command line arguments
p <- add_argument(p, "--input_path", help="Input path to seqtab object", type = "character")
p <- add_argument(p, "--output_path", help="Output path to final seqtab and taxo objects", type = "character")

# Parse the command line arguments
argv <- parse_args(p)
# print(argv)

# ----- READ DATA -----

seqtab.filename <- paste(argv$input_path, "/seqtab.rds", sep = "", collapse = TRUE)
seqtab <- readRDS(seqtab.filename)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
sum(seqtab.nochim)/sum(seqtab)

# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, "/users/work/cat3/db/dada2/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)

# Write to disk
ifelse(!dir.exists(file.path(argv$output_path)), dir.create(file.path(argv$output_path)), FALSE)
saveRDS(seqtab.nochim, paste(argv$output_path, "seqtab_final.rds", sep = "", collpase = TRUE)) # CHANGE ME to where you want sequence table saved
saveRDS(tax, paste(argv$output_path, "tax_final.rds", sep = "", collpase = TRUE))
