#!/usr/bin/env Rscript
# scripts/dada2-merge-chimera-taxo.r

# ./scripts/dada2-merge-chimera-taxo.r -i results/20190508_0074/dada2-inference/ -o results/20190508_0074/dada2-merge-chimera-taxo/

# ----- LIBRARIES -----

library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")

# library(argparser, quietly=TRUE)
#
# # ----- ARGUMENT PARSING -----
#
# # Create a parser
# p <- arg_parser("Run DADA2 Merge, remove chimera and assign taxo")
#
# # Add command line arguments
# p <- add_argument(p, "--input_path", help="Input path to seqtab object", type = "character")
# p <- add_argument(p, "--output_path", help="Output path to final seqtab and taxo objects", type = "character")
#
# # Parse the command line arguments
# argv <- parse_args(p)
# # print(argv)

# ----- READ DATA -----

seqtab.filename <- "/users/home/cat3/projects/mime-16s/results/20190503_0073/dada2-merge-chimera-taxo/seqtab_final.rds"

seqtab <- readRDS(seqtab.filename)

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))

fitGTR$tree
detach("package:phangorn", unload=TRUE)
