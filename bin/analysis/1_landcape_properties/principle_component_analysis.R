#!/usr/bin/env Rscript
# Perform principle component analysis on the full dataset

source('src/config.R')
source('src/dimensionality_reduction.R')

dir.create('figures/1_landscape_properties')
plots <- list()

#### Import Data ####
dms <- read_tsv('data/combined_mutational_scans.tsv')

foldx_averages <- select(dms, study, position, wt, total_energy:entropy_complex) %>%
  select(-sloop_entropy, -mloop_entropy, -entropy_complex, -water_bridge) %>% # Drop terms that are unused in our structures
  drop_na(total_energy) %>%
  group_by(study, position, wt) %>%
  summarise_all(mean, na.rm=TRUE)

position_constants <- select(dms, study, position, wt, phi:hydrophobicity) %>%
  distinct()

dms_wide <- filter(dms, mut %in% Biostrings::AA_STANDARD) %>%
  select(study, gene, position, wt, mut, imputed_score, log10_sift) %>%
  pivot_wider(names_from = mut, values_from = c(imputed_score, log10_sift)) %>%
  rename_at(vars(starts_with('imputed_score_')), ~str_sub(., start=-1))

pca <- tibble_pca(dms_wide, A:Y)
  
dms_wide <- bind_cols(dms_wide, as_tibble(pca$x)) %>%
  mutate(mean_score = rowMeans(select(., A:Y)),
         mean_sift = rowMeans(select(., log10_sift_A:log10_sift_Y))) %>%
  left_join(foldx_averages, by = c('study', 'position', 'wt')) %>%
  left_join(position_constants, by = c('study', 'position', 'wt'))
########

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

# FoldX term correlation
foldx_term_cor <- tibble_correlation(dms_wide, x = PC1:PC20, y = total_energy:energy_ionisation, filter_diag = TRUE, use = 'pairwise')
plots$foldx_pc_cor <- (ggplot(foldx_term_cor, aes(x=cat1, y=cat2, fill=cor)) +
  geom_raster() +
  scale_fill_gradient2(guide = guide_colourbar(title = 'Pearson\nCorrelation')) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))) %>%
  labeled_plot(units = 'cm', height = 15, width = 20)

# Save plots
save_plotlist(plots, 'figures/1_dimensionality_reduction/', overwrite = 'all', default_format = 'pdf')
