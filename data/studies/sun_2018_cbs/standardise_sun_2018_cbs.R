#!/usr/bin/env Rscript
# Standardise data from Sun et al. 2018 (CBS) (Preprint)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/sun_2018_cbs/sun_2018_cbs.yaml')
dm_data <- read_mavedb('data/studies/sun_2018_cbs/raw/urn_mavedb_00000005-a-4_scores.csv', score_transform = function(x){log2(x + min(x[x > 0], na.rm = TRUE))})

# Save output
standardise_study(dm_data, meta$study, meta$transform)