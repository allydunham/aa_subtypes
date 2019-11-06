#!/usr/bin/env Rscript
# Standardise data from Hartman et al. 2018 (CP)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/hartman_2018_cp/hartman_2018_cp.yaml')
dm_data <- read_xlsx('data/studies/hartman_2018_cp/raw/41467_2018_3783_MOESM4_ESM.xlsx', skip = 1, na = 'Not Observed') %>%
  rename(position = `Residue #`) %>%
  pivot_longer(-position, names_to = 'mut', values_to = 'raw_score') %>%
  drop_na(raw_score) %>% # Not all variants measured
  mutate(position = position + 1,
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut),
         raw_score = ifelse(raw_score == -4, -2.5, raw_score), # Set nulls to a value closer to rest of scale (they already set this arbitarily)
         transformed_score = raw_score/log10(2), # Transform base
         transformed_score = transformed_score - mean(transformed_score[class == 'Synonymous'], na.rm=TRUE), # 'divide' (in log domain) by WT scores
         score = normalise_score(transformed_score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
