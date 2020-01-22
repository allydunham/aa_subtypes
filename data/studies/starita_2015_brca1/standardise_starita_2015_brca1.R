#!/usr/bin/env Rscript
# Standardise data from Starita et al. 2015 (BRCA1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/starita_2015_brca1/starita_2015_brca1.yaml')
dm_data <- read_xls('data/studies/starita_2015_brca1/raw/genetics.115.175802-6.xls', na = 'NA') %>%
  rename_all(tolower) %>%
  rename(position = pos) %>%
  # Ref seq given by study has a mysterious, undocumented R at pos 175 where normal refs have K
  # using K here since the change is not explained in the paper and appears erroneous
  mutate(wt = str_split(meta$seq, '')[[1]][position], 
         class = get_variant_class(wt, mut)) %>%
  filter(!variant_id == 'NA-NA') %>%
  # Use E3 score - this is more general funcition based and empirically it seems to find most of the same negative effects as the BARD1 binding assay
  mutate(raw_score = e3_score, 
         transformed_score = transform_vamp_seq(raw_score),
         score = normalise_score(transformed_score)) %>%
  drop_na(score) # Not all measured

# Save output
standardise_study(dm_data, meta$study, meta$transform)
