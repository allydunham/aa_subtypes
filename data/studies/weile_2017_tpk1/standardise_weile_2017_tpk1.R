#!/usr/bin/env Rscript
# Standardise data from Weile et al. 2017 (TPK1)
source('src/config.R')
source('src/study_standardising.R')

aa_code <- structure(names(Biostrings::AMINO_ACID_CODE), names=Biostrings::AMINO_ACID_CODE)

# Import and process data
meta <- read_yaml('data/studies/weile_2017_tpk1/weile_2017_tpk1.yaml')
dm_data <- read_csv('data/studies/weile_2017_tpk1/raw/weile_2017_tpk1_score_comp.csv', na = c('NA','','None')) %>%
  mutate(mut = str_replace_all(hgvs_pro, aa_code)) %>%
  tidyr::extract(mut, into = c('wt', 'position', 'mut'), "p\\.([A-Z]+)([0-9]*)([A-Z=]+)", convert=TRUE) %>%
  mutate(mut = if_else(mut == '=', wt, mut),
         raw_score = score,
         score = transform_vamp_seq(raw_score),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, raw_score, class) %>%
  arrange(position, mut)

# Save output
standardise_study(dm_data, meta$study, meta$transform)