#!/usr/bin/env Rscript
# Script to validate the workflow used for Giacomelli et al. 2018 (TP53)
source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data_properties/per_study/giacomelli_2018_tp53')

# Import data
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
  mutate(fitness = log2(er / mean(er[variant_group == 'BackboneWt'], na.rm = TRUE)),
         class = get_variant_class(wt, mut)) %>%
  ungroup()

# Check whether the two experiments correlate
p_exp <- select(dm_data, -initial_freq, -freq, -er) %>%
  pivot_wider(names_from = experiment, values_from = fitness, names_prefix = 'Experiment ') %>%
  ggplot(aes(x = `Experiment 1`, y = `Experiment 2`)) +
  facet_grid(rows = vars(p53), cols = vars(drug)) +
  geom_point(aes(colour = class), shape = 20) +
  geom_smooth(colour = 'black', method = 'lm') +
  geom_abline(slope = 1, linetype = 'dashed') +
  labs(title = 'Correlation between technical replicates', subtitle = 'Without codon averaging')
ggsave('figures/0_data_properties/per_study/giacomelli_2018_tp53/initial_experiment_cor.pdf', p_exp, width = 15, height = 15, units = 'cm')

# Averaging codons improves cor between experiments
p_exp_codon <- group_by(dm_data, position, wt, mut, p53, experiment, drug, class) %>%
  summarise(fitness = weighted.mean(fitness, initial_freq, na.rm = TRUE)) %>% # weight by initial_freq to give more weight to measurements we're more confident of
  pivot_wider(names_from = experiment, values_from = fitness, names_prefix = 'Experiment ') %>%
  ggplot(aes(x = `Experiment 1`, y = `Experiment 2`)) +
  facet_grid(rows = vars(p53), cols = vars(drug)) +
  geom_point(aes(colour = class), shape = 20) +
  geom_smooth(colour = 'black', method = 'lm') +
  geom_abline(slope = 1, linetype = 'dashed') +
  labs(title = 'Correlation between technical replicates', subtitle = 'With codon averaging')
ggsave('figures/0_data_properties/per_study/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf', p_exp_codon, width = 15, height = 15, units = 'cm')

## Average codons and then experiments
dm_data <- group_by(dm_data, position, wt, mut, class, p53, experiment, drug) %>%
  summarise(raw_score = weighted.mean(fitness, initial_freq, na.rm = TRUE)) %>%
  group_by(position, wt, mut, class, p53, drug) %>%
  summarise(raw_score = mean(raw_score, na.rm = TRUE))

# Test three conditions:
p_cond <- ggplot(dm_data, aes(x = raw_score, fill = class)) +
  facet_grid(rows = vars(p53), cols = vars(drug)) +
  geom_histogram() +
  labs(title = 'Distribution of the three conditions in Giacomelli et al. 2018 (TP53)')
ggsave('figures/0_data_properties/per_study/giacomelli_2018_tp53/conditions.pdf', p_cond, width = 15, height = 15, units = 'cm')

# Etoposide both matches our desired distribution and selection type (enriches for functional TP53) so just use that
# others do test other aspects of the protein though and could potentially be integrated


