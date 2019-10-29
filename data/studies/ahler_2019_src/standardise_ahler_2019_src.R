#!/usr/bin/env Rscript
# Standardise data from Ahler et al. 2019 (Src)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/ahler_2019_src/ahler_2019_src.yaml')
dm_data <- bind_rows(read_mavedb('data/studies/ahler_2019_src/raw/urn_mavedb_00000041-b-1_scores.csv', position_offset = 1, score_col = activity_score),
                     read_mavedb('data/studies/ahler_2019_src/raw/urn_mavedb_00000041-a-1_scores.csv', position_offset = 269, score_col = activity_Score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)