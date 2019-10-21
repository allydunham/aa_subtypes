#!/usr/bin/env Rscript
# Standardise data from Firnberg et al. 2014 (TEM1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/firnberg_2014_tem1/firnberg_2014_tem1.yaml')
dm_data <- read_xlsx('data/studies/firnberg_2014_tem1/raw/firnberg_2014_tem1.xlsx') %>%
  rename(position = `Ambler Position`,
         ref_codon = `WT codon`,
         alt_codon = `Mutant codon`,
         wt = `WT AA`,
         mut = `Mutant AA`,
         base_changes = `Base Changes`,
         seq_counts_0.25 = `Sequencing Counts`,
         seq_counts_0.5 = ...8,
         seq_counts_1 = ...9,
         seq_counts_2 = ...10,
         seq_counts_4 = ...11,
         seq_counts_8 = ...12,
         seq_counts_16 = ...13,
         seq_counts_32 = ...14,
         seq_counts_64 = ...15,
         seq_counts_128 = ...16,
         seq_counts_256 = ...17,
         seq_counts_512 = ...18,
         seq_counts_1024 = ...19,
         total_seq_count = `Total Counts`,
         raw_score = Fitness,
         fitness_err = `Estimated error in fitness`) %>%
  drop_na(position) %>%
  filter(!wt == '*') %>%
  mutate(position = rep(1:nchar(meta$seq), each=64)) %>% # Numbering seems broken - starts at 3 and then misses 237 & 251
  group_by(position, wt, mut) %>%
  summarise(raw_score = mean(raw_score, na.rm = TRUE),
            score = mean(log2(score), na.rm = TRUE)) %>% # Average over codons
  mutate(score = normalise_score(score), 
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
