#!/usr/bin/env Rscript
# Standardise data from Dorrity et al. 2018 (STE12)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/dorrity_2018_ste12/dorrity_2018_ste12.yaml')
mating_data <- read_xlsx('data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd01.xlsx') %>%
  mutate(mut = sapply(seqID, process_split_seqid, USE.NAMES = FALSE),
         nmut = str_count(mut, ',') + 1,
         mating_avg = (mating_30C_rep1 + mating_30C_rep2 + mating_30C_rep3)/3) %>%
  select(mut, nmut, mating_avg, starts_with('mating_30C_'))

invasion_data <- read_xlsx('data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd02.xlsx') %>%
  mutate(mut = sapply(seqID, process_split_seqid, USE.NAMES = FALSE),
         nmut = str_count(mut, ',') + 1,
         invasion_avg = (invasion_30C_rep1 + invasion_30C_rep2 + invasion_30C_rep3)/3) %>%
  select(mut, nmut, invasion_avg, starts_with('invasion_30C_'))

dm_data <- full_join(mating_data, invasion_data, by = c('mut', 'nmut')) %>%
  
  # Take worst of two scores (both are essential functions)
  mutate(raw_score = pmin(mating_avg, invasion_avg, na.rm = TRUE)) %>%
  select(mut, nmut, raw_score) %>%
  separate(mut, into = str_c('mut', 1:max(.$nmut)), sep=',', fill = 'right') %>%
  pivot_longer(starts_with('mut'), names_to = 'tmp', values_to = 'mut') %>%
  select(-tmp) %>%
  drop_na(mut) %>%
  tidyr::extract(mut, into = c('position', 'mut'), '([0-9]*)([A-Z*])', convert = TRUE) %>%
  mutate(position = position + 140,
         wt = str_split(meta$seq, '')[[1]][position]) %>%

  # Take value from a single variant if possible, average over multiple mutations otherwise (<=4 gives 97% coverage without using very heavily mutated seqs)
  group_by(position, wt, mut) %>%
  summarise(raw_score = ifelse(any(nmut == 1), mean(raw_score[nmut == 1], na.rm = TRUE), mean(raw_score[nmut <= 1], na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(score = normalise_score(raw_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
