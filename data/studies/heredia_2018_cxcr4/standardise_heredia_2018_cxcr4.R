#!/usr/bin/env Rscript
# Standardise data from Heredia et al. 2018 (CXCR4)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/heredia_2018_cxcr4/heredia_2018_cxcr4.yaml')
dm_data <- read_xlsx('data/studies/heredia_2018_cxcr4/raw/GSE100368_enrichment_ratios_CXCR4.xlsx', skip = 6,
                     col_names = c('wt', 'mut', 'reads',
                                   'surface_exp_fitc_r1', 'surface_exp_fitc_r2',
                                   'binding_12g5_r1', 'binding_12g5_r2',
                                   'surface_exp_alexa_r1', 'surface_exp_alexa_r2',
                                   'binding_cxcl12_r1', 'binding_cxcl12_r2')) %>%
  mutate(wt = rep(wt[!is.na(wt)], each = 21)) %>%
  tidyr::extract(wt, into = c('wt', 'position'), '([A-Z])([0-9]+)', convert = TRUE) %>%
  
  # Average replicates for binding (which also incorporate surface expression)
  mutate(binding_12g5 = rowMeans(select(., binding_12g5_r1, binding_12g5_r2), na.rm = TRUE) %>% replace_na(NA),
         binding_cxcl12 = rowMeans(select(., binding_cxcl12_r1, binding_cxcl12_r2), na.rm = TRUE) %>% replace_na(NA)) %>%
  
  # Average binding scores for two conditions (0.48 correlation)
  mutate(raw_score = rowMeans(select(., binding_12g5, binding_cxcl12), na.rm = TRUE) %>% replace_na(NA),
         transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, transformed_score, raw_score, class) %>%
  drop_na(score) # not all measured

# Save output
standardise_study(dm_data, meta$study, meta$transform)
