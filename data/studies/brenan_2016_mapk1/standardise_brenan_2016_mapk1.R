#!/usr/bin/env Rscript
# Standardise data from Brenan et al. 2016 (MAPK1)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/brenan_2016_mapk1/brenan_2016_mapk1.yaml')
dm_data <- read_xlsx('data/studies/brenan_2016_mapk1/raw/brenan_2016_erk2.xlsx', sheet = 'Supplemental_Table_1') %>%
  rename_all(list( ~ gsub(' ', '_', tolower(.)))) %>%
  rename(wt = wt_aa, mut = mutant_aa, position = erk2_residue) %>%
  mutate(raw_score = `lfc_(etp_vs._dox)`, # Only general condition, other two are for specific drugs
         score = normalise_score(-raw_score), # The selection scheme they used favoured lof > wt > gof
         class = get_variant_class(wt, mut)) %>% 
  mutate_at(vars(nuc_acid_changes, dox_rank, sch_rank, vrt_rank, vrt_specific_allele, sch_specific_allele), as.integer) %>%
  select(position, wt, mut, score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)
