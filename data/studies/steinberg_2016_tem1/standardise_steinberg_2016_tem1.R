#!/usr/bin/env Rscript
# Standardise data from Steinberg & Ostermeier 2016 (TEM1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/steinberg_2016_tem1/steinberg_2016_tem1.yaml')
dm_data <- read_xlsx('data/studies/steinberg_2016_tem1/raw/1-s2.0-S0022283616301450-mmc2.xlsx', trim_ws = TRUE) %>%
  rename_all(~str_replace_all(str_to_lower(.), ' ', '_')) %>%
  drop_na(codon_position) %>%
  select(position=codon_position, wt = wt_aa, mut=mutant_aa, raw_score=tem1_amp_fitness) %>%
  mutate(score = normalise_score(log2(raw_score)),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)