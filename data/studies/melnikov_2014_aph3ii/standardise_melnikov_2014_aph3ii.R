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
meta <- read_yaml('IMPORT META YAML')
dm_data <- read_csv('IMPORT STUDY DATA')

# Save output
standardise_study(dm_data, meta$study, meta$transform)