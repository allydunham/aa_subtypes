#!/usr/bin/env Rscript
# Validate multi mutant combination method for Araya et al. 2012 (YAP1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/araya_2012_yap1/araya_2012_yap1.yaml')
dm_data <- read_csv('data/studies/araya_2012_yap1/raw/urn_mavedb_00000002-a-2_scores.csv', skip = 4) %>%
  select(hgvs_pro, raw_score = score) %>%
  mutate(hgvs_pro = if_else(str_ends(hgvs_pro, ']'), str_sub(hgvs_pro, start = 4, end = -2), str_sub(hgvs_pro, start = 3)),
         n_mut = str_count(hgvs_pro, ';') + 1) %>%
  separate(hgvs_pro, str_c('mut', 1:max(.$n_mut)), sep = ';', fill = 'right') %>%
  pivot_longer(cols = starts_with('mut'), values_to = 'mut') %>%
  drop_na(mut) %>%
  select(-name) %>%
  tidyr::extract(mut, into = c('wt', 'position', 'mut'), "([A-Za-z]{3})([0-9]+)([A-Za-z]{3})", convert = TRUE) %>%
  mutate(wt = AA_THREE_2_ONE[wt], mut = AA_THREE_2_ONE[mut], position = position + 169) %>%
  mutate(transformed_score = raw_score) %>%
  group_by(position, wt, mut)

singles <- summarise(dm_data, single_score := ifelse(any(n_mut == 1), mean(transformed_score[n_mut == 1], na.rm = TRUE), NA)) %>%
  ungroup()
single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)

mut_count_data <- lapply(2:9, function(x){
  summarise(dm_data, !!str_c('score_', x) := ifelse(any(n_mut <= x), mean(transformed_score[n_mut <= x], na.rm = TRUE), NA)) %>%
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
  geom_text(x = -5, y = 1, aes(label = str_c('frac = ', signif(frac, digits = 4))), hjust = 0) +
  labs(x = 'log2(Score) (Single Variant)',
       y = 'log2(Score) (Mean Over Multiple Variants)',
       title = 'Accuracy of multi-variant averaging for scoring in Araya et al. 2012 (YAP1)',
       subtitle = str_c('Fraction of variants with individual measures: ', single_frac))
ggsave('figures/0_data/per_study/araya_2012_yap1/multi_mut_validation.pdf', p_sing_multi, units = 'cm', height = 15, width = 25)
