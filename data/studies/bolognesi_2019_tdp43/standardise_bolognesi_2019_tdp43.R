#!/usr/bin/env Rscript
# Standardise data from Bolognesi et al. 2019 (TDP43)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/bolognesi_2019_tdp43/bolognesi_2019_tdp43.yaml')
dm_data <- read_xlsx('data/studies/bolognesi_2019_tdp43/raw/41467_2019_12101_MOESM7_ESM.xlsx', sheet = '1 AA change') %>%
  rename(position = Pos_abs, wt = WT_AA, mut = Mut, raw_score = toxicity) %>%
  mutate(transformed_score = raw_score / log(2),
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, transformed_score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
