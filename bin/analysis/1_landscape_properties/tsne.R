#!/usr/bin/env Rscript
# Perform tSNE analysis on the combined landscape

source('src/config.R')
source('src/dimensionality_reduction.R')
library(Rtsne)

dir.create('figures/1_landscape_properties')
plots <- list()
dms_wide <- read_tsv('data/combined_mutational_scans.tsv')

#### Analyse ####
plots$tsne_study <- (ggplot(dms_wide, aes(x=tSNE1, y=tSNE2, colour=study)) +
  geom_point() +
  facet_wrap(~study, labeller = as_labeller(sapply(unique(dms_wide$study), format_study, max_width=30))) +
  guides(colour=FALSE)) %>%
  labeled_plot(units='cm', height = 25, width = 25)

plots$tsne_aa <- (mutate(dms_wide, aa_class = AA_REDUCED_HASH[wt]) %>%
  ggplot(aes(x=tSNE1, y=tSNE2, colour=wt)) +
  geom_point() +
  facet_wrap(~aa_class) +
  scale_colour_manual(values = AA_COLOURS) +
  guides(colour = guide_legend(title='AA')) +
  theme(panel.grid.major.x = element_line(colour = 'gray', linetype = 'dotted'))) %>%
  labeled_plot(units='cm', height = 25, width = 25)

plots$tsne_hydrophobicity <- (ggplot(dms_wide, aes(x=tSNE1, y=tSNE2, colour=hydrophobicity)) +
  geom_point() +
  scale_colour_gradient2() +
  guides(colour = guide_colourbar(title = 'Hydrophobicity'))) %>%
  labeled_plot(units='cm', height = 10, width = 15)
  
plots$tsne_surface_accessibility <- (drop_na(dms_wide, all_atom_abs) %>%
                                     ggplot(aes(x=tSNE1, y=tSNE2, colour=all_atom_abs)) +
                                       geom_point() +
                                       scale_colour_viridis_c() +
                                       guides(colour = guide_colourbar(title = 'Surface Accessibility\n(All Atom Abs)'))) %>%
  labeled_plot(units='cm', height = 10, width = 15)


########

# Save plots
save_plotlist(plots, 'figures/1_landscape_properties/', overwrite = 'all', default_format = 'pdf')
