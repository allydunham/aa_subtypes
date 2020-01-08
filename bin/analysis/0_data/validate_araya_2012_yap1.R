#!/usr/bin/env Rscript
# Validate multi mutant combination method for Araya et al. 2012 (YAP1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/araya_2012_yap1/araya_2012_yap1.yaml')
dm_data <- read_tsv('data/studies/araya_2012_yap1/raw/araya_2012_hYAP65_ww.tsv', na = 'na',
                    col_types = cols(positions=col_character())) %>%
  separate(mutations, str_c('mutant', 1:max(.$mutation.count)), sep = ',', fill = 'right') %>%
  rename(raw_score = fitness, n_mut = mutation.count) %>%
  mutate(transformed_score = log2(raw_score)) %>%
  pivot_longer(starts_with('mutant'), values_to = 'mutant') %>%
  select(-name) %>%
  drop_na(mutant) %>%
  separate(mutant, into = c('position', 'mut'), sep = -1, convert = TRUE) %>%
  mutate(position = position + 160,
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut)) %>%
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
