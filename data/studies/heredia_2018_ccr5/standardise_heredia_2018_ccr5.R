#!/usr/bin/env Rscript
# Standardise data from Heredia et al. 2018 (CCR5)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/heredia_2018_ccr5/heredia_2018_ccr5.yaml')
dm_data <- read_xlsx('data/studies/heredia_2018_ccr5/raw/GSE100368_enrichment_ratios_CCR5.xlsx', skip = 7,
                     col_names = c('wt', 'mut', 'reads_l1', 'surface_exp_fitc_l1_r1', 'surface_exp_fitc_l1_r2', 'binding_2d7_l1_r1', 'binding_2d7_l1_r2',
                                   'empty', 'reads_l2', 'surface_exp_alexa_l2_r1', 'surface_exp_alexa_l2_r2', 'binding_gp120_cd4_l2_r1', 'binding_gp120_cd4_l2_r2')) %>%
  select(-empty) %>%
  mutate(wt = rep(wt[!is.na(wt)], each = 21)) %>%
  tidyr::extract(wt, into = c('wt', 'position'), '([A-Z])([0-9]+)', convert = TRUE) %>%
  
  # Average replicates for binding (which also incorporate surface expression)
  mutate(binding_2d7_l1 = rowMeans(select(., binding_2d7_l1_r1, binding_2d7_l1_r2), na.rm = TRUE) %>% replace_na(NA),
         binding_gp120_cd4_l2 = rowMeans(select(., binding_gp120_cd4_l2_r1, binding_gp120_cd4_l2_r2), na.rm = TRUE) %>% replace_na(NA)) %>%
  
  # Average binding scores for two conditions (0.61 correlation)
  mutate(raw_score = rowMeans(select(., binding_2d7_l1, binding_gp120_cd4_l2), na.rm = TRUE) %>% replace_na(NA),
         transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, transformed_score, raw_score, class) %>%
  drop_na(score) # not all measured


# Save output
standardise_study(dm_data, meta$study, meta$transform)
