#!/usr/bin/env Rscript
# Validate multiple mutation method for Starita et al. 2013 (UBE4B)
source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data/per_study/starita_2013_ube4b/')

# Import and process data
meta <- read_yaml('data/studies/starita_2013_ube4b/starita_2013_ube4b.yaml')
dm_data <- read_xlsx('data/studies/starita_2013_ube4b/raw/starita_2013_ube4b_ubox.xlsx', na = c('NA', '')) %>%
  filter(!seqID == 'NA-NA') %>% # Filter WT
  rename(raw_score = log2_ratio) %>%
  separate(seqID, into = c('position', 'mut'), sep='-') %>%
  select(-nscor_log2_ratio) %>%
  mutate(n_mut = sapply(position, function(x){str_count(x, ',') + 1})) %>%
  separate(mut, str_c('mut', 1:max(.$n_mut)), sep = ',', fill = 'right') %>%
  separate(position, str_c('position', 1:max(.$n_mut)), sep = ',', fill = 'right') %>%
  pivot_longer(starts_with('position'), names_to = 'pos_num', names_prefix = 'position', values_to = 'position') %>%
  drop_na(position) %>%
  pivot_longer(starts_with('mut'), names_to = 'mut_num', names_prefix = 'mut', values_to = 'mut') %>%
  drop_na(mut) %>%
  filter(pos_num == mut_num) %>%
  select(-pos_num, -mut_num) %>%
  group_by(position, mut)

singles <- summarise(dm_data, single_score := ifelse(any(n_mut == 1), mean(raw_score[n_mut == 1], na.rm = TRUE), NA)) %>%
  ungroup()
single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)

mut_count_data <- lapply(2:9, function(x){
  summarise(dm_data, !!str_c('score_', x) := ifelse(any(n_mut <= x), mean(raw_score[n_mut <= x], na.rm = TRUE), NA)) %>%
    ungroup() %>%
    select(!!str_c('score_', x))
}) %>%
  bind_cols(singles, .) %>%
  pivot_longer(starts_with('score_'), names_to = 'n_mut', values_to = 'multi_score', names_prefix = 'score_', names_ptypes = list(n_mut=integer())) %>%
  group_by(n_mut) %>%
  mutate(frac = sum(!is.na(multi_score))/length(multi_score)) %>%
  ungroup() %>%
  mutate(n_mut = as.integer(n_mut))

p_sing_multi <- ggplot(mut_count_data, aes(x = single_score, y = multi_score)) +
  facet_wrap(~n_mut, ncol = 5) +
  geom_point(colour = 'cornflowerblue') + 
  geom_abline(slope = 1, linetype='dashed') +
  geom_text(x = -4, y = 3, aes(label = str_c('frac = ', signif(frac, digits = 4))), hjust = 0) +
  labs(x = 'log2(Score) (Single Variant)',
       y = 'log2(Score) (Mean Over Multiple Variants)',
       title = 'Accuracy of multi-variant averaging for scoring in Starita et al. 2013 (UBE4B)',
       subtitle = str_c('Fraction of variants with individual measures: ', single_frac))
ggsave('figures/0_data/per_study/starita_2013_ube4b/multi_mut_validation.pdf', p_sing_multi, units = 'cm', height = 15, width = 25)
