#!/usr/bin/env Rscript
# Standardise data from Melnikov et al. 2014 (APH(3')-II)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/melnikov_2014_aph3ii/melnikov_2014_aph3ii.yaml')

count_files <- grep('\\.aacounts\\.txt', dir('data/studies/melnikov_2014_aph3ii/raw/'), value = TRUE)
count_files <- count_files[!grepl('(S[12]\\_Ami|S3\\_Kan)', count_files)] # Discard bad runs - see 00README.txt from Melnikov et al. data

counts <- sapply(count_files, read_melnikov_table, simplify = FALSE) %>%
  set_names(gsub('(KKA2\\_|\\.aacounts\\.txt)', '', names(.))) # Set names to drug
bkg_counts <- counts[c('Bkg1', 'Bkg2')]
counts <- counts[which(!names(counts) %in% c('Bkg1', 'Bkg2'))]

# Process data, see bin/0_data_properties/validate_melnikov.R for details + comments
dm_data <- mapply(melnikov_fitness, counts, names(counts), MoreArgs = list(bkg=bkg_counts), SIMPLIFY = FALSE) %>% # Calculate fitness for each count
  bind_rows(.id = 'experiment') %>%
  
  # Extract expeciment data, round and library contain the same information (plus round notes which needed a re-test at different MIC, 
  # which we alread accounted for) -> discard round
  separate(experiment, c('round', 'drug', 'library'), sep='_') %>%
  select(position, wt, mut, score, drug, library) %>% 
  
  # Process library pairs - discard datasets where libraries don't agree and filter outlier points, then average L1 & L2
  pivot_wider(id_cols = c('position', 'wt', 'mut', 'drug'), names_from = library, values_from = 'score') %>%
  mutate(diff = abs(L1 - L2)) %>%
  filter(dm_data, !(drug == 'Ami' & rel_conc == 0.25), !(drug %in% c('G418', 'Ami', 'Kan') & rel_conc == 0.125)) %>%
  filter(diff < sd(diff, na.rm = TRUE) * 3) %>%
  mutate(score = (L1 + L2)/2) %>%
  drop_na(score) %>%
  
  # Select the best conc for each drug, based on ER distribution
  filter(dm_data,
         (drug == 'Ami' & rel_conc == 0.5) |
           (drug == 'G418' & rel_conc == 0.25) |
           (drug == 'Kan' & rel_conc == 0.25) |
           (drug == 'Neo' & rel_conc == 0.25) |
           (drug == 'Paro' & rel_conc == 0.125) |
           (drug == 'Ribo' & rel_conc == 0.125)) %>%
  select(drug, position, wt, mut, score) %>%
  
  # Filter Ami as it doesn't correlate with other drugs, then average
  filter(dm_data, !drug == 'Ami') %>%
  group_by(position, wt, mut) %>%
  summarise(score = mean(score, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(raw_score = score,
         score = normalise_score(raw_score),
         class = get_variant_class(wt, mut))

# Save output
standardise_study(dm_data, meta$study, meta$transform)