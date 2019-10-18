#!/usr/bin/env Rscript
# Standardise data from Starita et al. 2015 (BRCA1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/starita_2015_brca1/starita_2015_brca1.yaml')
dm_data <- read_xls('data/studies/starita_2015_brca1/raw/starita_2015_brca1_ring.xls', na = 'NA') %>%
  rename_all(tolower) %>%
  rename(position = pos) %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position], # Ref seq given by study has a mysterious, undocumented R at pos 175 where normal refs have K, using K here since the change is not explained in the paper and appears erroneous
         class = get_variant_class(wt, mut)) %>%
  filter(!variant_id == 'NA-NA') %>%
  mutate(raw_score = pmin(e3_score, y2h_score, na.rm = TRUE),
         score = transform_vamp_seq(raw_score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)