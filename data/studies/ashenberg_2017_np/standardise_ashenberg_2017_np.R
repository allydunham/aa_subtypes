#!/usr/bin/env Rscript
# Standardise data from Ashenberg et al. 2017 (Flu NP)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/ashenberg_2017_np/ashenberg_2017_np.yaml')
dm_data <- read_csv('data/studies/ashenberg_2017_np/raw/journal.ppat.1006288.s013.csv') %>%
  rename(position = site, raw_score = diffsel) %>%
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  drop_na(score) # some muts just aren't measured

# Save output
standardise_study(dm_data, meta$study, meta$transform)
