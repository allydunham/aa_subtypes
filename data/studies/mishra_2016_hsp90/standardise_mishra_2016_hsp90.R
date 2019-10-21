#!/usr/bin/env Rscript
# Standardise data from Mishra et al. 2016 (HSP90)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/mishra_2016_hsp90/mishra_2016_hsp90.yaml')
path <- 'data/studies/mishra_2016_hsp90/raw/mishra_2016_hsp90_enrichment.xlsx'
dm_data <- map(excel_sheets(path), read_mishra_sheet, path = path) %>%
  bind_rows() %>%
  select(position, mut=aa, raw_score=avg) %>%
  mutate(raw_score = na_if(raw_score, -999),
         wt = str_split(meta$seq, '')[[1]][position],
         score = normalise_score(raw_score), 
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, raw_score, score, class) %>%
  arrange(position, mut)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
