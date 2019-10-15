#!/usr/bin/env Rscript
# Standardise data from Matreyek et al. 2018 (TPMT)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/metreyek_2018_tpmt/matreyek_2018_tpmt.yaml')
dm_data <- read_csv('data/studies/matreyek_2018_tpmt/raw/TPMT.csv',
                    col_types = cols(.default = col_character(), position = col_integer(), score = col_double())) %>%
  select(-X1) %>%
  rename(wt = start, mut = end, raw_score = score) %>%
  mutate(mut = if_else(mut == 'X', '*', mut),
         class = str_to_title(class),
         score = transform_vamp_seq(raw_score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

