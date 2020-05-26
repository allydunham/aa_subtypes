#!/usr/bin/env Rscript
# Produce figure S8-27 (Subtype Characterisation)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)

figures <- group_by(dms, wt) %>%
  group_map(~plot_full_characterisation(unique(.$cluster), full_characterisation, exclude_outliers = TRUE, global_scale = FALSE)) %>%
  map(extract2, 'overall') %>%
  set_names(str_c('figureS', 8:27))

save_plotlist(figures, 'figures/4_figures/', default_format = 'pdf')
save_plotlist(figures, 'figures/4_figures/', default_format = 'png')
