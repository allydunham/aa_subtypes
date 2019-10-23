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
  
  # Average replicates
  mutate(surface_exp_fitc_l1 = (surface_exp_fitc_l1_r1 + surface_exp_fitc_l1_r2)/2,
         binding_2d7_l1 = (binding_2d7_l1_r1 + binding_2d7_l1_r2)/2,
         surface_exp_alexa_l2 = (surface_exp_alexa_l2_r1 + surface_exp_alexa_l2_r2)/2,
         binding_gp120_cd4_l2 = (binding_gp120_cd4_l2_r1 + binding_gp120_cd4_l2_r2)/2,) %>%
  
  # Average surface expression score and take worst of three measured phenotypes
  mutate(surface_exp = (surface_exp_fitc_l1 + surface_exp_alexa_l2)/2,
         raw_score = pmin(binding_2d7_l1, surface_exp, binding_gp120_cd4_l2),
         score = normalise_score(raw_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, raw_score, class)


# Save output
standardise_study(dm_data, meta$study, meta$transform)