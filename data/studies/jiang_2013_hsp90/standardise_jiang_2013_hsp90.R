#!/usr/bin/env Rscript
# Standardise data from Jiang et al. 2013 (HSP90)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/jiang_2013_hsp90/jiang_2013_hsp90.yaml')
dm_data <- read_xlsx('data/studies/jiang_2013_hsp90/raw/jiang_2013_hsp90.xlsx', skip = 2) %>%
  select(-...10) %>%
  rename_all(tolower) %>%
  rename(mut = `amino acid`,
         sd = `standard deviation`) %>%
  mutate(raw_score = as.numeric(str_remove(average, '<')), # Set all instances of <0.034 to 0.034
         score = log2(raw_score),
         score = score / -min(score, na.rm = TRUE),
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut))
  
# Save output
standardise_study(dm_data, meta$study, meta$transform)

