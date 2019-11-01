#!/usr/bin/env Rscript
# Standardise data from Melamed et al. 2013 (PAB1)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/melamed_2013_pab1/melamed_2013_pab1.yaml')
dm_data <- read_xlsx('data/studies/melamed_2013_pab1/raw/melamed_2013_pab1_rrm_enrichment_ratios.xlsx') %>%
  rename(wt = WT_aa) %>%
  gather(key = 'mut', value = 'raw_score', -position, -wt) %>%
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

