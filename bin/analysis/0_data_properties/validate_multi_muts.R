#!/usr/bin/env Rscript
# Script testing performance of averaging over all occurances of a variant in multi datasets

source('src/config.R')

araya <- read_tsv('data/studies/araya_2012_yap1/raw/araya_2012_hYAP65_ww.tsv', na = 'na',
                             col_types = cols(positions=col_character())) %>%
  separate(mutations, str_c('mutant', 1:max(.$mutation.count)), sep = ',', fill = 'right') %>%
  rename(score = fitness, n = mutation.count, position = positions, mut = amino.acids) %>%
  mutate(gene = 'YAP1') %>%
  select(gene, position, mut, n, score)

starita <- read_xlsx('data/studies/starita_2013_ube4b/raw/starita_2013_ube4b_ubox.xlsx', na = c('NA', '')) %>%
  rename(score = log2_ratio) %>%
  separate(seqID, into = c('position', 'mut'), sep='-') %>%
  mutate(gene = 'UBE4B', n = sapply(position, function(x){str_count(x, ',') + 1})) %>%
  select(-nscor_log2_ratio)

variants <- bind_rows(araya, starita) %>%
  separate(mut, str_c('mut', 1:max(.$n)), sep = ',', fill = 'right') %>%
  separate(position, str_c('position', 1:max(.$n)), sep = ',', fill = 'right') %>%
  pivot_longer(starts_with('position'), names_to = 'pos_num', names_prefix = 'position', values_to = 'position') %>%
  drop_na(position) %>%
  pivot_longer(starts_with('mut'), names_to = 'mut_num', names_prefix = 'mut', values_to = 'mut') %>%
  drop_na(mut) %>%
  filter(pos_num == mut_num) %>%
  select(-pos_num, -mut_num) %>%
  group_by(gene, position, mut) %>%
  summarise(Mean_avg_score = mean(score, na.rm = TRUE),
            Median_avg_score = median(score, na.rm = TRUE),
            sing_score = ifelse(any(n == 1), score[which(n == 1)[1]], NA)) %>%
  pivot_longer(ends_with('_avg_score'), names_to = 'avg', names_pattern = "(.*)_avg_score", values_to = 'score')

p <- ggplot(variants, aes(x = sing_score, y = score, colour = gene)) +
  facet_wrap(~avg, ncol = 1) +
  geom_point() +
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = 'Score from Single Variants', y = 'Score from Average of Multiples',
       title = 'Assesment of imputing single variant scores by averaging multi-mutants',
       subtitle = 'Using Araya et al. 2012 (YAP1) & Starita et al. 2013 (UBE4B)') +
  guides(colour = guide_legend(title = '', override.aes = list(fill=NA)))

ggsave('figures/0_data_properties/averaging_multi_mutants.pdf', p, units = 'cm', height = 22, width = 18)
