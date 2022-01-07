library(microbiome) # Load libraries
library(phyloseq)
library(dplyr)
library(reshape2)
library(knitr)
library(devtools)
library(phylosmith)

# ----- PRELIMINARY SETTINGS -----

# import colours to remain consistent between all plots
source_url("https://raw.githubusercontent.com/clarajegousse/mime-16s/main/scripts/visual-settings.r")

# ----- LOAD PHYLOSEQ OBJ -----

ps <- readRDS("/Users/Clara/Projects/mime-16s/global-ps-ark.rds")

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# ----- SELECT PROKARYOTES ONLY -----

# because these were assigned with Silva
ps0 <- subset_taxa(ps, Kingdom %in% c("Archaea"))

# ----- PREVALENCE FILTERING -----

# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(ps0),
              MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps0)
ps1

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))
ps2

ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum) + theme_pubr() +
  #scale_color_manual(values = phylum.color) +
  tax_color_scale(ps0, "Phylum") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ----- AGGLOMERATE CLOSELY RELATED TAXA -----

# Taxonomic agglomeration
# How many genera are present after filtering?
taxGlomRank = "Order"
length(get_taxa_unique(ps2, taxonomic.rank = taxGlomRank))

ps3 = tax_glom(ps2, taxrank = taxGlomRank)
ps3

# ------- SUBSET SPECIFIC SAMPLES ----

# here let's select only the surface samples

#ps4 <- subset_samples(ps2, depth =="0")
ps4 <- ps3

# ----- SUBSET SPECIFIC GROUP -----

ps5 <- subset_taxa(ps4, Kingdom %in% c("Archaea"))
ps5

# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps5))
standf = function(x, t=total) round(t * (x / sum(x)))
#ps6 = transform_sample_counts(ps5, standf)
ps6 = transform_sample_counts(ps5, log10)

sum(is.na(otu_table(ps6)))
psmelt(ps6) %>%
  filter(is.na(Abundance))


abundance_heatmap(ps5, 
                  classification = 'Phylum', 
                  treatment = c('transect', 'zone'), 
                  subset = 'year',
                  )

