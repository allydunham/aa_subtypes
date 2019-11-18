#!/usr/bin/env Rscript
# Combine all standardised datasets

source('src/config.R')
source('src/study_standardising.R')

study_dirs <- commandArgs(trailingOnly = TRUE)

dms <- sapply(study_dirs, import_study, fields = c('gene'), simplify = FALSE) %>%
  bind_rows()

# Import Sift results
sift <- sapply(unique(dms$gene), import_sift, simplify = FALSE) %>%
  bind_rows(.id = 'gene') 

# Import FoldX results
foldx <- sapply(unique(dms$gene), import_foldx, simplify = FALSE) %>%
  bind_rows(.id = 'gene')

# Combine data
dms <- left_join(dms, sift, by = c('gene', 'position', 'wt', 'mut')) %>%
  left_join(foldx, by = c('gene', 'position', 'wt', 'mut'))

# Write combined data
write_tsv(dms, 'data/combined_mutational_scans.tsv')
