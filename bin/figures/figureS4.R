#!/usr/bin/env Rscript
# Produce figure S4 (Study Confounding)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv') %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))

umap2_breaks <- c(-2.5, 0, 2.5)


