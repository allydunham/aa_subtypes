#!/usr/bin/env Rscript
# Perform principle component analysis on the full dataset

source('src/config.R')
source('src/dimensionality_reduction.R')

dir.create('figures/1_landscape_properties')
plots <- list()
dms_wide <- read_tsv('data/combined_mutational_scans.tsv')

#### Plot PCs vs factors ####
# PC1 vs significance
plots$pc1_vs_mean_score <- (ggplot(dms_wide, aes(x = mean_score, y = PC1, colour = mean_sift)) +
  geom_point() +
  labs(x = 'Mean Normalised ER', y = 'PC1') +
  scale_color_gradient(low = 'black', high = 'red', guide = guide_colourbar(title = 'Mean log10(SIFT)'))) %>%
  labeled_plot(units = 'cm', height = 8, width = 12)

plots$pc1_pc2_mean_score <- (ggplot(dms_wide, aes(x = PC1, y = PC2, colour = mean_score)) +
                              geom_point() +
                              labs(x = 'PC1', y = 'PC2') +
                              scale_color_gradient(low = 'black', high = 'red', guide = guide_colourbar(title = 'Mean Normalised ER'))) %>%
  labeled_plot(units = 'cm', height = 8, width = 12)

# PCs vs surface accessibility
plots$pc2_pc4_surface_accessibility <- (ggplot(dms_wide, aes(x = PC2, y = PC4, colour = all_atom_abs)) +
                                        geom_point() +
                                        labs(x = 'PC2', y = 'PC4') +
                                        scale_color_viridis_c(guide = guide_colourbar(title = 'Surface Accessibility'))) %>%
  labeled_plot(units = 'cm', height = 8, width = 12)

plots$pc2_vs_surface_accessibility <- (ggplot(dms_wide, aes(x = PC2, y = polar_abs, colour = hydrophobicity)) +
                                         geom_point() +
                                         labs(x = 'PC2', y = 'Polar Residue Surface Accesibility') +
                                         scale_color_viridis_c(guide = guide_colourbar(title = 'Hydrophobicity'))) %>%
  labeled_plot(units = 'cm', height = 8, width = 12)

plots$pc2_vs_hydrophobicity <- (ggplot(dms_wide, aes(x = hydrophobicity, y = PC2)) +
                                  geom_boxplot(aes(group=cut(hydrophobicity, 10))) +
                                  geom_smooth(method = 'lm') +
                                  labs(y = 'PC2', x = 'Hydrophobicity')) %>%
  labeled_plot(units = 'cm', height = 8, width = 12)

# PC3 vs FoldX terms
plots$pc3_foldx_terms <- (ggplot(dms_wide, aes(x = entropy_sidechain, y=van_der_waals, colour=PC3)) + 
                            geom_point() + 
                            scale_colour_distiller(type = 'div', palette = 'Spectral', values = c(0, 0.4, 0.5, 0.6, 1)) + 
                            labs(x = expression('Sidechain Entropy (kcal mol'^-1*')'), y = expression('Van der Waals (kcal mol'^-1*')'))) %>%
  labeled_plot(units = 'cm', height = 10, width = 15)

# FoldX term correlation
foldx_term_cor <- tibble_correlation(dms_wide, x = PC1:PC20, y = total_energy:energy_ionisation, filter_diag = TRUE, use = 'pairwise') %>%
  mutate(cat1 = factor(cat1, levels = levels(cat1)[order(as.integer(str_sub(levels(cat1), start = 3)))]))
plots$foldx_pc_cor <- (ggplot(foldx_term_cor, aes(x=cat1, y=cat2, fill=cor)) +
  geom_raster() +
  scale_fill_gradient2(guide = guide_colourbar(title = 'Pearson\nCorrelation')) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))) %>%
  labeled_plot(units = 'cm', height = 15, width = 20)

# Save plots
save_plotlist(plots, 'figures/1_landscape_properties/', overwrite = 'all', default_format = 'pdf')
