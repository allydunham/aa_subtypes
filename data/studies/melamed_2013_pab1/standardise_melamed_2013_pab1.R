#!/usr/bin/env Rscript
# Standardise data from Melamed et al. 2013 (PAB1)

source('src/config.R')
source('src/study_standardising.R')

study_id = 'melamed_2013_pab1'
transform = 'None'

# Import and process data
meta <- read_yaml('data/studies/melamed_2013_pab1/melamed_2013_pab1.yaml')
dm_data <- read_xlsx('data/studies/melamed_2013_pab1/raw/melamed_2013_pab1_rrm_enrichment_ratios.xlsx') %>%
  rename(wt = WT_aa) %>%
  gather(key = 'mut', value = 'raw_score', -position, -wt) %>%
  mutate(score = raw_score / -min(raw_score, na.rm = TRUE), class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, study_id, transform)