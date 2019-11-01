#!/usr/bin/env Rscript
# Standardise data from Araya et al. 2012 (YAP1)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/araya_2012_yap1/araya_2012_yap1.yaml')
dm_data <- read_tsv('data/studies/araya_2012_yap1/raw/araya_2012_hYAP65_ww.tsv', na = 'na',
                    col_types = cols(positions=col_character())) %>%
  separate(mutations, str_c('mutant', 1:max(.$mutation.count)), sep = ',', fill = 'right') %>%
  rename(raw_score = fitness) %>%
  mutate(transformed_score = log2(raw_score)) %>%
  pivot_longer(starts_with('mutant'), values_to = 'mutant') %>%
  select(-name) %>%
  drop_na(mutant) %>%
  separate(mutant, into = c('position', 'mut'), sep = -1, convert = TRUE) %>%
  group_by(position, mut) %>%
  summarise(transformed_score = mean(transformed_score), raw_score = mean(raw_score)) %>% # Average over occurances of a variant?
  ungroup() %>%
  mutate(position = position + 160,
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut),
         score = normalise_score(score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
