#!/usr/bin/env Rscript
# Standardise data from Doud et al. 2015 (NP)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/doud_2015_np/doud_2015_np.yaml')
dm_data <- read_table2('data/studies/doud_2015_np/raw/Supp_file_2_mean_aichi1968_prefs.txt', skip = 1,
                       col_names = c('position', 'wt', 'entropy', 'A', 'C', 'D', 'E',
                                     'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                     'R', 'S', 'T', 'V', 'W', 'Y'),
                       col_types = cols(position=col_integer(), wt=col_character(), .default = col_double())) %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'raw_score') %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position]) %>% # is the same in all but 334, where the data has '?', so replace with known WTs
  group_by(position) %>%
  mutate(transformed_score = log2(raw_score / raw_score[which(mut == first(wt))])) %>% # Normalise by the WT at that position
  ungroup() %>%
  mutate(score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
