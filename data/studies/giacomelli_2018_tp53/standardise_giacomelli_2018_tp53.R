#!/usr/bin/env Rscript
# Standardise data from Giacomelli et al. (TP53)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/giacomelli_2018_tp53/giacomelli_2018_tp53.yaml')
dm_data <- read_xlsx('data/studies/giacomelli_2018_tp53/raw/41588_2018_204_MOESM5_ESM.xlsx', skip=1) %>%
  rename_all(str_to_lower) %>%
  rename(wt = wt_aa, mut = vt_aa) %>%
  mutate(mut = if_else(mut == 'Z', '*', mut),
         wt = if_else(wt == 'Z', '*', wt)) %>%
  filter(a549_p53wt_early_time_point_experiment_1 > 0,
         a549_p53wt_early_time_point_experiment_2 > 0) %>%
  
  # Transform to frequencies with pseudocount
  mutate_at(vars(starts_with('a549')), .funs = ~(. + min(.[. > 0], na.rm = TRUE))/sum(., na.rm = TRUE)) %>%
  
  # Extract columns for the combination of states
  pivot_longer(starts_with('a549'), names_to = 'id', values_to = 'count') %>%
  tidyr::extract(id, into = c('p53', 'state', 'experiment'), "a549_p53(wt|null)_(early_time_point|nutlin-3|etoposide)_experiment_([12])") %>%
  pivot_wider(names_from = state, values_from = count) %>%
  rename(nutlin3 = `nutlin-3`, initial_freq = early_time_point) %>%
  pivot_longer(c('nutlin3', 'etoposide'), names_to = 'drug', values_to = 'freq') %>%
  drop_na(freq) %>%
  
  # Determine ER and fitness for each combination
  mutate(er = freq / initial_freq) %>%
  select(position, wt, mut, wt_codon, vt_codon, variant_group, p53, experiment, drug, initial_freq, freq, er) %>%
  group_by(p53, experiment, drug) %>%
  mutate(fitness = log2(er / mean(er[variant_group == 'BackboneWt'], na.rm = TRUE))) %>%
  
  # Average codons and then experiments
  group_by(position, wt, mut, p53, experiment, drug) %>%
  summarise(raw_score = weighted.mean(fitness, initial_freq, na.rm = TRUE)) %>%
  group_by(position, wt, mut, p53, drug) %>%
  summarise(raw_score = mean(raw_score, na.rm = TRUE)) %>%
  ungroup() %>%
  
  # Select appropriate experiment (p53 NULL, Etoposide selects for funcional p53), others test other functions that could possibly be integrated carefully
  filter(p53 == 'null', drug == 'etoposide') %>%
  
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score), 
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, transformed_score, raw_score, class)
  

# Save output
standardise_study(dm_data, meta$study, meta$transform)
