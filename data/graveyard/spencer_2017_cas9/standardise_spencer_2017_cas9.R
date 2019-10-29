#!/usr/bin/env Rscript
# Standardise data from Spencer et al. 2017 (Cas9)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/spencer_2017_cas9/spencer_2017_cas9.yaml')
dm_data <- read_xlsx('data/studies/spencer_2017_cas9/raw/41598_2017_17081_MOESM2_ESM.xlsx') %>%
  rename_all(~str_replace_all(str_to_lower(.), ' ', '_')) %>%
  select(-...37) %>%
  select(position = aa_position, wt = wt_aa, mut= mutant_aa, raw_score = log2_fold_change_after_positive_selection) %>%
  mutate(score = normalise_score(raw_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)