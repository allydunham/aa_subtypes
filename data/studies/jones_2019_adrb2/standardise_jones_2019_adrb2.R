#!/usr/bin/env Rscript
# Standardise data from Jones et al. 2019 (ADRB2)

source('src/config.R')
source('src/study_standardising.R')

study_id = 'jones_2019_adrb2'
transform = 'log2(x/2)'

# Import and process data
meta <- read_yaml('data/studies/jones_2019_adrb2/jones_2019_adrb2.yaml')
dm_data <- read_csv('data/studies/jones_2019_adrb2/raw/lib-med.csv') %>%
  filter(Condition == 0.150) %>% # Select only EC50 measure (generally correlate between these)
  select(position = Pos, mut = AA, raw_score = Median, Repeat) %>%
  group_by(position, mut) %>%
  summarise(raw_score = median(raw_score)) %>% # Average biological repeats
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         score = log2(raw_score / 2), # Divide by 2 as twice missense peaks around 2?
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, study_id, transform)