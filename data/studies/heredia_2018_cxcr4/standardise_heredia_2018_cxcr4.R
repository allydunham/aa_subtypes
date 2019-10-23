#!/usr/bin/env Rscript
# Standardise data from STUDY

## README - delete in real script versions
# Template script for standardising study data
# Saves a standardised tsv /data/studies/{study}/{study}.tsv
# With columns: position, wt, mut, score (which is the normalised score, all this processing is done here)
# Transform score such that -1 = NULL, 0 = WT, +ve = beneficial
# Tibble passed to standardise_study() must have at least columns position, wt, mut, score, raw_score

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
  
  # Average replicates
  mutate(surface_exp_fitc = (surface_exp_fitc_r1 + surface_exp_fitc_r2)/2,
         binding_12g5 = (binding_12g5_r1 + binding_12g5_r2)/2,
         surface_exp_alexa = (surface_exp_alexa_r1 + surface_exp_alexa_r2)/2,
         binding_cxcl12 = (binding_cxcl12_r1 + binding_cxcl12_r2)/2,) %>%
  
  # Average surface expression score and take worst of three measured phenotypes
  mutate(surface_exp = (surface_exp_fitc + surface_exp_alexa)/2,
         raw_score = pmin(binding_12g5, surface_exp, binding_cxcl12),
         score = normalise_score(raw_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)