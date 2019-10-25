#!/usr/bin/env Rscript
# Standardise data from Adkar et al. 2012 (CcdB)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/adkar_2012_ccdb/adkar_2012_ccdb.yaml')
dm_data <- read_xls('data/studies/adkar_2012_ccdb/raw/1-s2.0-S0969212612000068-mmc2.xls', skip = 1) %>%
  rename_all(~str_to_lower(str_replace(., ' ', '')))

# Save output
standardise_study(dm_data, meta$study, meta$transform)