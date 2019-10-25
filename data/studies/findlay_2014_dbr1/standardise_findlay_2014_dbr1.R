#!/usr/bin/env Rscript
# Standardise data from Findlay et al. 2014 (DBR1
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/findlay_2014_dbr1/findlay_2014_dbr1.yaml')
dm_data <- read_xlsx('data/studies/findlay_2014_dbr1/raw/findlay_2014_dbr1_exon2_counts.xlsx', skip = 3, na = 'NA') %>%
  rename(seq = Sequence,
         log2_enrichment_score_day11_rep1 = `Day 11 log2 enrichment score (replicate 1)`,
         log2_enrichment_score_day11_rep2 = `Day 11 log2 enrichment score (replicate 2)`,
         position = `Affected codon`,
         wt = `Reference AA`,
         mut = `Substituted AA`,
         mut_type = `Variant Class`) %>%
  mutate(position = as.integer(position),
         mut = ifelse(mut == 'WT', wt, mut),
         raw_score = (log2_enrichment_score_day11_rep1 + log2_enrichment_score_day11_rep2)/2,
         score = normalise_score(raw_score),
         class = get_variant_class(wt, mut)) %>%
  drop_na(position, raw_score) %>%
  filter(!mut == 'DEL') %>%
  select(position, wt, mut, score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)