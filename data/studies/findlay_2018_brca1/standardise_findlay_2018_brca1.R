#!/usr/bin/env Rscript
# Standardise data from Findlay et al. 2018 (BRCA1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/findlay_2018_brca1/findlay_2018_brca1.yaml')
dm_data <- read_xlsx('data/studies/findlay_2018_brca1/raw/findlay_2018_brca1_ring_brct.xlsx', skip = 2, na = 'NA') %>%
  rename_all(list( ~ gsub('[\\/ \\(\\)]+', '_', .))) %>%
  rename(wt_nuc = reference,
         mut_nuc = alt,
         wt = aa_ref,
         mut = aa_alt,
         position = aa_pos) %>%
  drop_na(position) %>%
  mutate(raw_score = function.score.mean,
         transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, transformed_score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
