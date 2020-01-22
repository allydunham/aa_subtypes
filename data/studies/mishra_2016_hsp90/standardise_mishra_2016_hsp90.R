#!/usr/bin/env Rscript
# Standardise data from Mishra et al. 2016 (HSP90)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/mishra_2016_hsp90/mishra_2016_hsp90.yaml')
path <- 'data/studies/mishra_2016_hsp90/raw/1-s2.0-S2211124716303175-mmc2.xlsx'
dm_data <- map(excel_sheets(path), read_mishra_sheet, path = path) %>%
  bind_rows() %>%
  select(position, mut=aa, raw_score=avg) %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         raw_score = na_if(raw_score, -999),
         transformed_score = raw_score,
         score = normalise_score(transformed_score), 
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, transformed_score, raw_score, score, class) %>%
  arrange(position, mut) %>%
  drop_na(score) # Some not measured

# Save output
standardise_study(dm_data, meta$study, meta$transform)
