#!/usr/bin/env Rscript
# Validate approach to Heredia et al. 2018 (CCR5 and CXCR4)
source('src/config.R')
source('src/study_standardising.R')

#### CCR5 ####
dir.create('figures/0_data/per_study/heredia_2018_ccr5')
dm_data <- read_xlsx('data/studies/heredia_2018_ccr5/raw/GSE100368_enrichment_ratios_CCR5.xlsx', skip = 7,
                     col_names = c('wt', 'mut', 'reads_l1', 'surface_exp_fitc_l1_r1', 'surface_exp_fitc_l1_r2', 'binding_2d7_l1_r1', 'binding_2d7_l1_r2',
                                   'empty', 'reads_l2', 'surface_exp_alexa_l2_r1', 'surface_exp_alexa_l2_r2', 'binding_gp120_cd4_l2_r1', 'binding_gp120_cd4_l2_r2')) %>%
  select(-empty) %>%
  mutate(wt = rep(wt[!is.na(wt)], each = 21)) %>%
  tidyr::extract(wt, into = c('wt', 'position'), '([A-Z])([0-9]+)', convert = TRUE)

# Library Correlation
p_ccr5_lib_cor <- select(dm_data, -starts_with('reads_')) %>%
  pivot_longer(surface_exp_fitc_l1_r1:binding_gp120_cd4_l2_r2, names_to = 'exp', values_to = 'er') %>%
  tidyr::extract(exp, into = c('condition', 'library', 'replicate'), "([a-z0-9\\_]*)_l([12])_r([12])") %>%
  pivot_wider(names_from = replicate, values_from = er, names_prefix = 'er_rep_') %>%
  ggplot(aes(x = er_rep_1, y = er_rep_2)) +
  facet_wrap(~condition) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  geom_abline(slope = 1, linetype = 'dotted') +
  labs(title = 'Correlation between replicates in Heredia et al. 2018 (CCR5)',
       x = 'Rep 1 ER', y = 'Rep 2 ER')
ggsave('figures/0_data/per_study/heredia_2018_ccr5/replicate_correlation.pdf', p_ccr5_lib_cor, units = 'cm', width = 20, height = 20)

# Doesn't agree hugely well, but average...
dm_data <- mutate(dm_data,
                  surface_exp_fitc_l1 = (surface_exp_fitc_l1_r1 + surface_exp_fitc_l1_r2)/2,
                  binding_2d7_l1 = (binding_2d7_l1_r1 + binding_2d7_l1_r2)/2,
                  surface_exp_alexa_l2 = (surface_exp_alexa_l2_r1 + surface_exp_alexa_l2_r2)/2,
                  binding_gp120_cd4_l2 = (binding_gp120_cd4_l2_r1 + binding_gp120_cd4_l2_r2)/2,) %>%
  select(position, wt, mut, surface_exp_fitc_l1:binding_gp120_cd4_l2)

# Correlation between experiments
p_ccr5_exp_cor <- ggpairs(dm_data, columns = 4:7)
ggsave('figures/0_data/per_study/heredia_2018_ccr5/experiment_correlation.pdf', p_ccr5_exp_cor, units = 'cm', width = 20, height = 20)

# Average surface expression tests, then take worst as they measure different aspects of fitness
dm_data <- mutate(dm_data,
                  surface_exp = (surface_exp_fitc_l1 + surface_exp_alexa_l2)/2,
                  raw_score = pmin(binding_2d7_l1, surface_exp, binding_gp120_cd4_l2),
                  score = normalise_score(raw_score),
                  class = get_variant_class(wt, mut))
########

#### CXCR4 ####
dir.create('figures/0_data/per_study/heredia_2018_cxcr4')
dm_data <- read_xlsx('data/studies/heredia_2018_cxcr4/raw/GSE100368_enrichment_ratios_CXCR4.xlsx', skip = 6,
                     col_names = c('wt', 'mut', 'reads',
                                   'surface_exp_fitc_r1', 'surface_exp_fitc_r2',
                                   'binding_12g5_r1', 'binding_12g5_r2',
                                   'surface_exp_alexa_r1', 'surface_exp_alexa_r2',
                                   'binding_cxcl12_r1', 'binding_cxcl12_r2')) %>%
  mutate(wt = rep(wt[!is.na(wt)], each = 21)) %>%
  tidyr::extract(wt, into = c('wt', 'position'), '([A-Z])([0-9]+)', convert = TRUE)

# Library Correlation
p_cxcr4_lib_cor <- select(dm_data, -reads) %>%
  pivot_longer(surface_exp_fitc_r1:binding_cxcl12_r2, names_to = 'exp', values_to = 'er') %>%
  tidyr::extract(exp, into = c('condition', 'replicate'), "([a-z0-9\\_]*)_r([12])") %>%
  pivot_wider(names_from = replicate, values_from = er, names_prefix = 'er_rep_') %>%
  ggplot(aes(x = er_rep_1, y = er_rep_2)) +
  facet_wrap(~condition) +
  geom_point() +
  geom_smooth(method = 'lm') + 
  geom_abline(slope = 1, linetype = 'dotted') +
  labs(title = 'Correlation between replicates in Heredia et al. 2018 (CXCR4)',
       x = 'Rep 1 ER', y = 'Rep 2 ER')
ggsave('figures/0_data/per_study/heredia_2018_cxcr4/replicate_correlation.pdf', p_cxcr4_lib_cor, units = 'cm', width = 20, height = 20)

# Doesn't agree hugely well, but average...
dm_data <- mutate(dm_data,
                  surface_exp_fitc = (surface_exp_fitc_r1 + surface_exp_fitc_r2)/2,
                  binding_12g5 = (binding_12g5_r1 + binding_12g5_r2)/2,
                  surface_exp_alexa = (surface_exp_alexa_r1 + surface_exp_alexa_r2)/2,
                  binding_cxcl12 = (binding_cxcl12_r1 + binding_cxcl12_r2)/2,) %>%
  select(position, wt, mut, surface_exp_fitc:binding_cxcl12)

# Correlation between experiments
p_cxcr4_exp_cor <- ggpairs(dm_data, columns = 4:7)
ggsave('figures/0_data/per_study/heredia_2018_cxcr4/experiment_correlation.pdf', p_cxcr4_exp_cor, units = 'cm', width = 20, height = 20)

# Average surface expression tests, then take worst as they measure different aspects of fitness
dm_data <- mutate(dm_data,
                  surface_exp = (surface_exp_fitc + surface_exp_alexa)/2,
                  raw_score = pmin(binding_12g5, surface_exp, binding_cxcl12),
                  score = normalise_score(raw_score),
                  class = get_variant_class(wt, mut))
########