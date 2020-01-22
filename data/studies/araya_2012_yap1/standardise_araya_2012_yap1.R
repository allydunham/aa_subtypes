#!/usr/bin/env Rscript
# Standardise data from Araya et al. 2012 (YAP1)
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
  group_by(position, wt, mut) %>%
  summarise(transformed_score = ifelse(1 %in% n_mut, mean(transformed_score[n_mut == 1], na.rm=TRUE), mean(transformed_score[n_mut <= 2], na.rm=TRUE)),
            raw_score = ifelse(1 %in% n_mut, mean(raw_score[n_mut == 1], na.rm=TRUE), mean(raw_score[n_mut <= 2], na.rm=TRUE))) %>%
  ungroup() %>%
  mutate(class = get_variant_class(wt, mut),
         score = normalise_score(transformed_score)) %>%
  drop_na(score) %>% # Some mutant arent found in seqs with <= 2 variants
  arrange(position, mut)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
