#!/usr/bin/env Rscript
# Validate the strategy used for Sarkisyan et al. 2016 (GFP)
source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data_properties/per_study/sarkisyan_2016_gfp')

# Import and process data
raw_data <- read_tsv('data/studies/sarkisyan_2016_gfp/raw/sarkisyan_2016_gfp_AAs.tsv', skip = 1,
                     col_names = c('mut', 'barcodes', 'median_brightness', 'std')) %>%
            mutate(n_mut = str_count(mut, ':') + 1)
wt_brightness <- filter(raw_data, is.na(mut)) %>% pull(median_brightness)

dm_data <- separate(raw_data, mut, into = str_c('mut', 1:15), sep = ':', fill = 'right') %>%
    pivot_longer(cols = starts_with('mut'), names_to = 'n', names_prefix = 'mut', values_to = 'mut') %>%
    drop_na(mut) %>%
    select(-n, -barcodes, -std) %>%
    tidyr::extract(mut, into = c('wt', 'position', 'mut'), 'S([A-Z])([0-9]+)([A-Z])', convert=TRUE) %>%
    arrange(position, mut) %>%
    mutate(single_score = ifelse(n_mut == 1, median_brightness, NA)) %>%
    group_by(position, wt, mut)

singles <- summarise(dm_data, single_score := ifelse(any(n_mut == 1), mean(median_brightness[n_mut == 1], na.rm = TRUE), NA)) %>%
  ungroup()
single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)

mut_count_data <- lapply(2:15, function(x){
  summarise(dm_data, !!str_c('score_', x) := ifelse(any(n_mut <= x), mean(median_brightness[n_mut <= x], na.rm = TRUE), NA)) %>%
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
  geom_text(x = 1.5, y = 3.5, aes(label = str_c('frac = ', signif(frac, digits = 4))), hjust = 0) +
  labs(x = 'Median Brightness (Single Variant)',
       y = 'Median Brightness (Mean Over Multiple Variants)',
       title = 'Accuracy of multi-variant averaging for scoring in Sarkisyan et al. 2016 (GFP)',
       subtitle = str_c('Fraction of variants with individual measures: ', single_frac))
ggsave('figures/0_data_properties/per_study/sarkisyan_2016_gfp/multi_mut_validation.pdf', p_sing_multi, units = 'cm', height = 15, width = 25)
