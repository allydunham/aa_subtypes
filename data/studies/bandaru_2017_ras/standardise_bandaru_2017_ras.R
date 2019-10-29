#!/usr/bin/env Rscript
# Standardise data from Bandaru et al. 2017 (Ras)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/bandaru_2017_ras/bandaru_2017_ras.yaml')
dm_data <- read_xlsx('data/studies/bandaru_2017_ras/raw/elife-27810-supp1-v2.xlsx') %>%
  transpose_tibble(col_names = ...1, id_col = 'position')

# Save output
standardise_study(dm_data, meta$study, meta$transform)