#!/usr/bin/env Rscript
# Standardise data from Lee et al. 2018 (HA)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/lee_2018_ha/lee_2018_ha.yaml')

# Extract overall protein position - is split into signal peptide and HA1/HA2 chains in original file
get_position <- function(x){
  if (grepl('\\-', x)){
    return(as.numeric(x) + 17)
  } else if (grepl('HA2', x)) {
    x <- gsub('\\(HA2\\)', '', x)
    return(as.numeric(x) + 345)
  } else {
    return(as.numeric(x) + 16)
  }
}

dm_data <- read_xlsx('data/studies/lee_2018_ha/raw/pnas.1806133115.sd03.xlsx', sheet = 'avg_prefs') %>%
  rename(position = site) %>%
  gather(key = 'mut', value = 'raw_score', -position, -entropy, -neffective) %>%
  mutate(position = sapply(position, get_position),
         wt = str_split(meta$seq, '')[[1]][position],
         class = get_variant_class(wt, mut)) %>%
  arrange(position, mut) %>%
  group_by(position) %>%
  mutate(transformed_score = log2(raw_score / raw_score[which(mut == first(wt))])) %>%
  ungroup() %>%
  mutate(score = normalise_score(transformed_score))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
