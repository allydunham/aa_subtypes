#!/usr/bin/env Rscript
# Combine all standardised datasets

source('src/config.R')
source('src/study_standardising.R')

study_dirs <- commandArgs(trailingOnly = TRUE)

dms <- sapply(study_dirs, import_study, fields = c('gene'), simplify = FALSE) %>%
  bind_rows() %>%
  group_by(study, gene, position, wt) %>%
  filter(sum(!mut == wt) >= 15) %>% # Only keep positions with a maximum of 4 missing scores
  ungroup()

# Import Sift results
sift <- sapply(unique(dms$gene), import_sift, simplify = FALSE) %>%
  bind_rows(.id = 'gene') 

# Import FoldX results
foldx <- sapply(unique(dms$gene), import_foldx, simplify = FALSE) %>%
  bind_rows(.id = 'gene')

# Filter incomplete positions and impute
imputed_dms <- select(dms, -transformed_score, -raw_score) %>%
  pivot_wider(id_cols = c(study, gene, position, wt, mut), names_from = mut, values_from = score) %>%
  pivot_longer(c(-study, -position, -wt, -gene), names_to = 'mut', values_to = 'score') %>%
  group_by(wt, mut) %>%
  mutate(score = ifelse(is.na(score), ifelse(wt == mut, 0, median(score, na.rm = TRUE)), score)) %>%
  ungroup()

# Combine data
dms <- rename(imputed_dms, imputed_score=score) %>%
  left_join(dms, by = c('study', 'position', 'wt', 'mut', 'gene')) %>%
  left_join(sift, by = c('gene', 'position', 'wt', 'mut')) %>%
  left_join(foldx, by = c('gene', 'position', 'wt', 'mut')) %>%
  mutate(class = get_variant_class(wt, mut)) %>%
  arrange(study, position, wt, mut)

# Write combined data
write_tsv(dms, 'data/combined_mutational_scans.tsv')
