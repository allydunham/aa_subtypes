#!/usr/bin/env Rscript
# Standardise data from Kitzman et al. 2015 (GAL4)

source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/kitzman_2015_gal4/kitzman_2015_gal4.yaml')
path <- 'data/studies/kitzman_2015_gal4/raw/kitzman_2015_gal4_enrichment.xlsx'
dm_data <- lapply(excel_sheets(path), read_kitzman_sheet, path = path) %>%
  bind_rows(.) %>%
  spread(key = 'label', value = 'log2_enrichment') %>%
  mutate(raw_score = (SEL_A_24h + SEL_A_40h + SEL_B_40h + SEL_C_40h + SEL_C_64h)/5, # Average over replicates (diff times found to still correlate well)
         score = raw_score / -min(raw_score, na.rm = TRUE),
         class = get_variant_class(wt, mut)) %>%
  filter(!mut == 'delInFrame')

# Save output
standardise_study(dm_data, meta$study, meta$transform)