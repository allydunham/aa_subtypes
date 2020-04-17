#!/usr/bin/env Rscript
# Produce figure 4 (C Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

er_limits <- c(min(full_characterisation$profiles$er), -min(full_characterisation$profiles$er))

plot_profiles <- function(clusters, legend = FALSE){
  filter(full_characterisation$profiles, cluster %in% clusters) %>%
    mutate(mut = add_markdown(mut, AA_COLOURS),
           cluster = add_markdown(cluster, cluster_number_colourmap(cluster))) %>%
    ggplot(aes(x = mut, y = cluster, fill = er)) +
    geom_raster(show.legend = legend) +
    coord_fixed() +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
    guides(fill = guide_colourbar(title = 'Normalised ER')) + 
    theme(axis.text.x = element_markdown(),
          axis.text.y = element_markdown(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major.y = element_blank())
}

p_cys_profiles <- plot_profiles(c('C1', 'C2'))

p_cys_surface_acc <- filter(dms, cluster %in% c('C1', 'C2')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density', show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_cys_disulphide <- filter(dms, cluster %in% c('C1', 'C2'), !is.na(disulfide)) %>%
  ggplot(aes(x = cluster, fill = cluster)) +
  geom_bar(show.legend = FALSE) +
  labs(y = 'Disulphide Bonds', x = '') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  theme(legend.title = element_blank(),
        panel.grid.major.y = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))
er_legend <- plot_profiles('A1', legend = TRUE) %>% get_legend() %>% as_ggplot() + size

figure4 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size +
  labs(tag = 'A') +
  annotation_custom(ggplotGrob(p_cys_profiles + size), xmin = 0.15, xmax = 0.8, ymin = 0.4, ymax = 0.6) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_disulphide.png'), interpolate = TRUE,
                    xmin = 0.025, xmax = 0.2, ymin = 0.7, ymax = 0.875) +
  annotation_custom(ggplotGrob(p_cys_disulphide + size), xmin = 0.2, xmax = 0.475, ymin = 0.7, ymax = 0.875) +
  annotate('text', x = 0.25, y = 0.65, label = str_wrap('Most disulphide bonds are C1 positions', 25)) + 
  annotation_custom(ggplotGrob(p_cys_surface_acc + size), xmin = 0.525, xmax = 0.975, ymin = 0.7, ymax = 0.975) +
  annotate('text', x = 0.75, y = 0.65, label = str_wrap('C2 positions are more buried', 25)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/gal4_cys_zinc.png'), interpolate = TRUE,
                    xmin = 0.1, xmax = 0.4, ymin = 0.1, ymax = 0.4) +
  annotate('text', x = 0.25, y = 0.05, label = str_wrap('Other active roles are generally C1 positions', 25)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_aromatic.png'), interpolate = TRUE,
                    xmin = 0.525, xmax = 0.73, ymin = 0.15, ymax = 0.355) +
  annotation_raster(readPNG('figures/4_figures/position_examples/np_cys_aromatic.png'), interpolate = TRUE,
                    xmin = 0.77, xmax = 0.975, ymin = 0.15, ymax = 0.355) +
  annotate('text', x = 0.75, y = 0.05, label = str_wrap('Both C1 & C2 positions have aromatic interactions', 30))

ggsave('figures/4_figures/figure4.pdf', figure4, width = 200, height = 200, units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = 200, height = 200, units = 'mm')
