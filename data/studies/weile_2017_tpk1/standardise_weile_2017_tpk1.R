#!/usr/bin/env Rscript
# Standardise data from Weile et al. 2017 (TPK1)
source('src/config.R')
source('src/study_standardising.R')

aa_code <- structure(names(Biostrings::AMINO_ACID_CODE), names=Biostrings::AMINO_ACID_CODE)

# Import and process data
meta <- read_yaml('data/studies/weile_2017_tpk1/weile_2017_tpk1.yaml')
dm_data <- read_mavedb('data/studies/weile_2017_tpk1/raw/urn_mavedb_00000001-d-2_scores.csv', score_transform = transform_vamp_seq)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
