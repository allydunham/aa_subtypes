#!/usr/bin/env Rscript
# Standardise data from Mishra et al. 2016 (HSP90)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/mishra_2016_hsp90/mishra_2016_hsp90.yaml')
path <- 'data/studies/mishra_2016_hsp90/raw/mishra_2016_hsp90_enrichment.xlsx'
dm_data <- map(excel_sheets(path), read_mishra_sheet, path = path) %>%
  bind_rows() %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         score = raw_score / -min(raw_score, na.rm = TRUE),
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, raw_score, score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)