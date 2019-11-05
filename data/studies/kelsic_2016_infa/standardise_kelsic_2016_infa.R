#!/usr/bin/env Rscript
# Standardise data from Kelsic et al. 2016 (infA)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/kelsic_2016_infa/kelsic_2016_infa.yaml')
dm_data <- read_csv('data/studies/kelsic_2016_infa/raw/cels_206_mmc5.csv', skip = 1,
                    col_names = c('codon', 'mut', 'position', 'is_wt', 'raw_score',
                                  'sd', 'fitness_rich',
                                  'fitness_stdev_rich', 'RCS',
                                  'mfe_ddG_43nt_sliding', 'tmp')) %>%
  select(codon:sd) %>%
  mutate(transformed_score = transform_vamp_seq(raw_score)) %>%
  group_by(position, mut) %>%
  summarise_at(vars(transformed_score, raw_score), mean, na.rm=TRUE) %>%
  ungroup() %>%
  mutate(score = normalise_score(transformed_score),
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut)) %>%
  drop_na(wt, position) # measured codons for stop codon, which we don't want

# Save output
standardise_study(dm_data, meta$study, meta$transform)
