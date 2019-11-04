#!/usr/bin/env Rscript
# Standardise data from Starita et al. 2013 (UBE4B)

source('src/config.R')
source('src/study_standardising.R')

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
  group_by(position, mut) %>%
  summarise(raw_score = ifelse(1 %in% n_mut, mean(raw_score[n_mut == 1], na.rm=TRUE), mean(raw_score[n_mut <= 3], na.rm=TRUE))) %>%
  ungroup() %>%
  mutate(transformed_score = raw_score,
         score = normalise_score(transformed_score), 
         position = as.integer(position) + 1072, # tested region starts at +1072 according to Starita (slightly before uniprot UBOX) This does lead to ref seq aligning
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

