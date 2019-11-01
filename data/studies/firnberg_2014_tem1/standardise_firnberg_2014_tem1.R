#!/usr/bin/env Rscript
# Standardise data from Firnberg et al. 2014 (TEM1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/firnberg_2014_tem1/firnberg_2014_tem1.yaml')
dm_data <- read_xlsx('data/studies/firnberg_2014_tem1/raw/firnberg_2014_tem1.xlsx', skip = 1,
                     col_names = c('position', 'ref_codon', 'alt_codon', 'wt', 'mut', 'base_changes', 'seq_counts_0.25',
                                   'seq_counts_0.5', 'seq_counts_1', 'seq_counts_2', 'seq_counts_4', 'seq_counts_8',
                                   'seq_counts_16', 'seq_counts_32', 'seq_counts_64', 'seq_counts_128', 'seq_counts_256', 
                                   'seq_counts_512', 'seq_counts_1024', 'total_seq_count', 'raw_score', 'fitness_err')) %>%
  drop_na(position) %>%
  filter(!wt == '*') %>%
  mutate(position = rep(1:nchar(meta$seq), each=64)) %>% # Numbering seems broken - starts at 3 and then misses 237 & 251
  group_by(position, wt, mut) %>%
  summarise(raw_score = mean(raw_score, na.rm = TRUE),
            transformed_score = mean(log2(raw_score), na.rm = TRUE)) %>% # Average over codons
  mutate(score = normalise_score(transformed_score), 
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
