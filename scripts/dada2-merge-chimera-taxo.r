#!/usr/bin/env Rscript
# scripts/dada2-merge-chimera-taxo.r

# ./scripts/dada2-merge-chimera-taxo.r --input_path results/20190508_0074/dada2-inference/ --output_path results/20190508_0074/dada2-merge-chimera-taxo/

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

seqtab.filename <- paste(argv$input_path, "/seqtab.rds", sep = "", collapse = NULL)
seqtab <- readRDS(seqtab.filename)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
# sum(seqtab.nochim)/sum(seqtab)

# Choose database
db <- "GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz"
#db <- "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
db.path <- paste0("/users/work/cat3/db/dada2/", db)

# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, db.path, multithread=TRUE)

# Write to disk
ifelse(!dir.exists(file.path(argv$output_path)), dir.create(file.path(argv$output_path), showWarnings = FALSE), FALSE)

saveRDS(seqtab.nochim, paste(argv$output_path, "/seqtab_final.rds", sep = "", collapse = NULL)) # CHANGE ME to where you want sequence table saved
saveRDS(tax, paste(argv$output_path, "/tax_final.rds", sep = "", collapse = NULL))
