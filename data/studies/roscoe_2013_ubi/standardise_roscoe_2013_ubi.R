#!/usr/bin/env Rscript
# Standardise data from Roscoe et al. 2013 (Ubi) 

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/roscoe_2013_ubi/roscoe_2013_ubi.yaml')
dm_data <- read_xlsx('data/studies/roscoe_2013_ubi/raw/roscoe_2013_ubi_fitness.xlsx', skip = 4) %>%
  rename(position = Position,
         mut = `Amino Acid`,
         selection_chr = Apparent,
         sd_chr = `Quantified Synonyms`) %>%
  mutate(raw_score = as.numeric(selection_chr),
         score =normalise_score(raw_score), 
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, raw_score, score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)

