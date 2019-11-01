#!/usr/bin/env Rscript
# Standardise data from Bandaru et al. 2017 (Ras)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/bandaru_2017_ras/bandaru_2017_ras.yaml')
dm_data <- read_xlsx('data/studies/bandaru_2017_ras/raw/elife-27810-supp1-v2.xlsx') %>%
  mutate(...1 = replace_na(...1, 'wt')) %>%
  transpose_tibble(col_names = ...1, id_col = 'position') %>%
  mutate_at(vars(-wt), as.numeric) %>%
  mutate(position = as.integer(position)) %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'raw_score') %>%
  mutate(transformed_score = raw_score/log(2),
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
