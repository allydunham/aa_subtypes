#!/usr/bin/env Rscript
# Standardise data from Jones et al. 2019 (ADRB2)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/jones_2019_adrb2/jones_2019_adrb2.yaml')
dm_data <- read_csv('data/studies/jones_2019_adrb2/raw/lib-med.csv') %>%
  filter(Condition == 0.150) %>% # Select only EC50 measure (generally correlate between these)
  select(position = Pos, mut = AA, raw_score = Median, Repeat) %>%
  group_by(position, mut) %>%
  summarise(raw_score = median(raw_score)) %>% # Average biological repeats
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         transformed_score = log2(raw_score / 2), # Divide by 2 as neutral seems to peaks around 2?
         score = normalise_score(transformed_score), 
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

