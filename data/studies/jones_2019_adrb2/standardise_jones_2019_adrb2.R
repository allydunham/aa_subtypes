#!/usr/bin/env Rscript
# Standardise data from Jones et al. 2019 (ADRB2)

source('src/config.R')
source('src/study_standardising.R')

# Determine WT like
data("BLOSUM62")
similar_aas <- as_tibble(BLOSUM62, rownames = 'aa1') %>%
  pivot_longer(-aa1, names_to = 'aa2', values_to = 'blosum') %>%
  filter(blosum > 2, !aa1 == aa2, !aa1 %in% c('B', 'J', 'Z', 'X', '*'), !aa2 %in% c('B', 'J', 'Z', 'X', '*'))
similar_aas <- str_c(similar_aas$aa1, similar_aas$aa2)
  
# Import and process data
meta <- read_yaml('data/studies/jones_2019_adrb2/jones_2019_adrb2.yaml')
dm_data <- read_csv('data/studies/jones_2019_adrb2/raw/lib-med.csv') %>%
  filter(Condition == 0.150) %>% # Select only EC50 measure (generally correlate between these)
  select(position = Pos, mut = AA, raw_score = Median, Repeat) %>%
  group_by(position, mut) %>%
  summarise(raw_score = median(raw_score)) %>% # Average biological repeats
  ungroup() %>%
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         pair = str_c(wt, mut),
         # Divide by average score of v. similar AAs (blosum62 > 1) as substitute for wt
         transformed_score = log2(raw_score / mean(raw_score[pair %in% similar_aas], na.rm=TRUE)),
         score = normalise_score(transformed_score), 
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)

