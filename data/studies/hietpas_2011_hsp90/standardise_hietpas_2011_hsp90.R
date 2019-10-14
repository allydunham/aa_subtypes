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

study_id = 'heitpas_2011_hsp90'
transform = 'None'

# Import and process data
meta <- read_yaml('data/studies/hietpas_2011_hsp90/hietpas_2011_hsp90.yaml')

dm_data <- read_csv('data/studies/hietpas_2011_hsp90/raw/hietpas_2011_pdz_ligands_fitness.csv') %>%
  rename(mut = aa, raw_score = selection_coefficient) %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position]) %>%
  group_by(position, wt, mut) %>%
  summarise(raw_score = mean(raw_score)) %>%
  mutate(score = raw_score / -min(raw_score, na.rm = TRUE))

# Save output
standardise_study(dm_data, study_id, transform)