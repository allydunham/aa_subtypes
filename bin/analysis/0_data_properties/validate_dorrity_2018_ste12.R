#!/usr/bin/env Rscript
# Validate methodology used for Dorrity et al. 2018 (Ste12)

source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data_properties/per_study/dorrity_2018_ste12')

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

dm_data <- full_join(mating_data, invasion_data, by = c('mut', 'nmut'))

# Averaging reps
p_paired <- ggpairs(dm_data, columns = c('mating_avg', 'mating_30C_rep1', 'mating_30C_rep2', 'mating_30C_rep3',
                                         'invasion_avg', 'invasion_30C_rep1', 'invasion_30C_rep2', 'invasion_30C_rep3'))
ggsave('figures/0_data_properties/per_study/dorrity_2018_ste12/rep_correlation.pdf', p_paired, units = 'cm', height = 30, width = 30)

# Take worst of two traits, as both are protein functions
dm_data <- mutate(dm_data, raw_score = pmin(mating_avg, invasion_avg, na.rm = TRUE)) %>%
  select(mut, nmut, raw_score) %>%
  separate(mut, into = str_c('mut', 1:max(.$nmut)), sep=',', fill = 'right') %>%
  pivot_longer(starts_with('mut'), names_to = 'tmp', values_to = 'mut') %>%
  select(-tmp) %>%
  drop_na(mut) %>%
  tidyr::extract(mut, into = c('position', 'mut'), '([0-9]*)([A-Z*])', convert = TRUE) %>%
  mutate(position = position + 140,
         wt = str_split(meta$seq, '')[[1]][position])

singles <- group_by(dm_data, position, wt, mut) %>%
  summarise(single_score := ifelse(any(nmut == 1), mean(raw_score[nmut == 1], na.rm = TRUE), NA)) %>%
  ungroup()
single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)

mut_count_data <- lapply(2:max(dm_data$nmut), function(x){
  group_by(dm_data, position, wt, mut) %>%
  summarise(!!str_c('score_', x) := ifelse(any(nmut <= x), mean(raw_score[nmut <= x], na.rm = TRUE), NA)) %>%
    ungroup() %>%
    select(!!str_c('score_', x))
}) %>%
  bind_cols(singles, .) %>%
  pivot_longer(starts_with('score_'), names_to = 'nmut', values_to = 'multi_score', names_prefix = 'score_', names_ptypes = list(n_mut=integer())) %>%
  group_by(nmut) %>%
  mutate(frac = sum(!is.na(multi_score))/length(multi_score)) %>%
  ungroup() %>%
  mutate(nmut = as.integer(nmut))

p_sing_multi <- ggplot(mut_count_data, aes(x = single_score, y = multi_score)) +
  facet_wrap(~nmut, ncol = 5) +
  geom_point(colour = 'cornflowerblue') + 
  geom_abline(slope = 1, linetype='dashed') +
  geom_text(x = -6.25, y = 1.25, aes(label = str_c('frac = ', signif(frac, digits = 4))), hjust = 0) +
  labs(x = 'ER (Single Variant)',
       y = 'ER (Mean Over Multiple Variants)',
       title = 'Accuracy of multi-variant averaging for scoring in Dorrity et al. 2018 (STE12)',
       subtitle = str_c('Fraction of variants with individual measures: ', signif(single_frac, digits = 4)))
ggsave('figures/0_data_properties/per_study/dorrity_2018_ste12/multi_mut_validation.pdf', p_sing_multi, units = 'cm', height = 15, width = 25)
