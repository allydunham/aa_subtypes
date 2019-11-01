#!/usr/bin/env Rscript
# Standardise data from Wrenbeck et al. 2017 (amiE)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/wrenbeck_2017_amie/wrenbeck_2017_amie.yaml')
dm_data <- bind_rows(acetamide = read_tsv('data/studies/wrenbeck_2017_amie/raw/amiESelectionFitnessData_Acetamide.txt', na = c('NS', 'None')),
                     isobutyramide = read_tsv('data/studies/wrenbeck_2017_amie/raw/amiESelectionFitnessData_Isobutyramide.txt', na = c('NS', 'None')),
                     propionamide = read_tsv('data/studies/wrenbeck_2017_amie/raw/amiESelectionFitnessData_Propionamide.txt', na = c('NS', 'None')),
                     .id = 'condition') %>%
  select(position = location, mut = mutation, raw_score = normalized_fitness, condition) %>%
  pivot_wider(names_from = condition, values_from = raw_score) %>%
  mutate(raw_score = rowMeans(select(., acetamide, isobutyramide, propionamide), na.rm = TRUE) %>% replace_na(NA),
         transformed_score = raw_score,
         score = normalise_score(transformed_score),
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
