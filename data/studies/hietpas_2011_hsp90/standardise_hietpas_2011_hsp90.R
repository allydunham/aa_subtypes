#!/usr/bin/env Rscript
# Standardise data from Hietpas et al. 2011 (HSP90) 
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/hietpas_2011_hsp90/hietpas_2011_hsp90.yaml')
dm_data <- read_csv('data/studies/hietpas_2011_hsp90/raw/hietpas_2011_pdz_ligands_fitness.csv') %>%
  rename(mut = aa, raw_score = selection_coefficient) %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position]) %>%
  group_by(position, wt, mut) %>% # Average over codons
  summarise(raw_score = mean(raw_score)) %>%
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

