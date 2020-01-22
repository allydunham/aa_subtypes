#!/usr/bin/env Rscript
# Standardise data from Sarkisyan et al. 2016 (GFP)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/sarkisyan_2016_gfp/sarkisyan_2016_gfp.yaml')
raw_data <- read_tsv('data/studies/sarkisyan_2016_gfp/raw/amino_acid_genotypes_to_brightness.tsv', skip = 1,
                    col_names = c('mut', 'barcodes', 'median_brightness', 'std'))

wt_brightness <- filter(raw_data, is.na(mut)) %>% pull(median_brightness)

dm_data <- mutate(raw_data, n_mut = str_count(mut, ':') + 1) %>%
  filter(n_mut <= 3) %>%
  separate(mut, into = str_c('mut', 1:3), sep = ':', fill = 'right') %>%
  pivot_longer(cols = starts_with('mut'), names_to = 'n', names_prefix = 'mut', values_to = 'mut') %>%
  drop_na(mut) %>%
  select(-n, -barcodes, -std) %>%
  tidyr::extract(mut, into = c('wt', 'position', 'mut'), 'S([A-Z])([0-9]+)([A-Z*])', convert=TRUE, remove=FALSE) %>%
  mutate(position = position + 2) %>% # Numbered from 3rd residue for some reason
  arrange(position, mut) %>%
  group_by(position, wt, mut) %>%
  summarise(raw_score = if_else(1 %in% n_mut, mean(median_brightness[n_mut == 1], na.rm = TRUE), # Use value of single mut if available
                                mean(median_brightness[n_mut <= 4], na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(transformed_score = log2(raw_score / wt_brightness),
         score = normalise_score(transformed_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)
