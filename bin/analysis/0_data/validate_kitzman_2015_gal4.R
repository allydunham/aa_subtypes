#!/usr/bin/env Rscript
# Validate the selection experiment combination method for Kitzman et al. 2015 (GAL4)
source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data/per_study/kitzman_2015_gal4')

# Import and process data
path <- 'data/studies/kitzman_2015_gal4/raw/kitzman_2015_gal4_enrichment.xlsx'
dm_data <- lapply(excel_sheets(path), read_kitzman_sheet, path = path) %>%
  bind_rows(.) %>%
  spread(key = 'label', value = 'log2_enrichment') %>%
  filter(!mut == 'delInFrame')

# Plot variants
p <- ggpairs(dm_data, columns = c('NONSEL_24h', 'SEL_A_24h', 'SEL_A_40h', 'SEL_B_40h', 'SEL_C_40h', 'SEL_C_64h'))
ggsave('figures/0_data/per_study/kitzman_2015_gal4/validate_selection_combination.pdf', p, units = 'cm', width = 25, height = 25)