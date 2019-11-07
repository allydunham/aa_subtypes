#!/usr/bin/env Rscript
# Standardise data from Findlay et al. 2014 (DBR1
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/findlay_2014_dbr1/findlay_2014_dbr1.yaml')
dm_data <- read_xlsx('data/studies/findlay_2014_dbr1/raw/findlay_2014_dbr1_exon2_counts.xlsx', skip = 5, na = 'NA',
                     col_names = c('seq', 'log2_enrichment_score_day11_rep1', 'log2_enrichment_score_day11_rep2',
                                   'position', 'wt', 'mut', 'mut_type')) %>%
  mutate(position = as.integer(na_if(position, "3'SS")),
         mut = ifelse(mut == 'WT', wt, mut),
         raw_score = rowMeans(select(., log2_enrichment_score_day11_rep1, log2_enrichment_score_day11_rep2), na.rm = TRUE) %>% replace_na(NA)) %>%
  group_by(position, wt, mut) %>%
  summarise(raw_score = mean(raw_score, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut)) %>%
  drop_na(position, raw_score) %>%
  filter(!mut == 'DEL') %>%
  select(position, wt, mut, score, transformed_score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
