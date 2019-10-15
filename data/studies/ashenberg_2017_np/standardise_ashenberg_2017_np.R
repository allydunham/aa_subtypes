#!/usr/bin/env Rscript
# Standardise data from Ashenberg et al. 2017 (Flu NP)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_csv('data/studies/ashenberg_2017_np/ashenberg_2017_np.yaml')
dm_data <- read_csv('data/studies/ashenberg_2017_np/raw/ashenberg_2017_flu_np.csv') %>%
  rename(position = site, raw_score = diffsel) %>%
  mutate(score = raw_score / -min(raw_score, na.rm = TRUE),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
