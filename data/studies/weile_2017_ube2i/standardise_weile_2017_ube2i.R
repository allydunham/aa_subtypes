#!/usr/bin/env Rscript
# Standardise data from Weile et al. 2017 (UBE2I)
source('src/config.R')
source('src/study_standardising.R')

aa_code <- structure(names(Biostrings::AMINO_ACID_CODE), names=Biostrings::AMINO_ACID_CODE)

# Import and process data
meta <- read_yaml('data/studies/weile_2017_ube2i/weile_2017_ube2i.yaml')
dm_data <- read_mavedb('data/studies/weile_2017_ube2i/raw/urn_mavedb_00000001-a-1_scores.csv', score_transform = transform_vamp_seq) %>%
  filter(!position == 159) # Drop position not in WT protein

# Save output
standardise_study(dm_data, meta$study, meta$transform)
