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

study_id = 'roscoe_2013_ubi4'
transform = 'None'

# Import and process data
meta <- read_yaml('data/studies/roscoe_2013_ubi/roscoe_2013_ubi.yaml')
dm_data <- read_xlsx('data/studies/roscoe_2013_ubi/raw/roscoe_2013_ubi_fitness.xlsx', skip = 4) %>%
  rename(position = Position,
         mut = `Amino Acid`,
         selection_chr = Apparent,
         sd_chr = `Quantified Synonyms`) %>%
  mutate(raw_score = as.numeric(selection_chr),
         score = raw_score / -min(raw_score, na.rm = TRUE),
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, raw_score, score, class)

# Save output
standardise_study(dm_data, study_id, transform)