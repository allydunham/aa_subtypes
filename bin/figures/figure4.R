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

p_profiles <- filter(full_characterisation$profiles, cluster %in% c('C1', 'C2')) %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         cluster = add_markdown(cluster, cluster_number_colourmap(cluster))) %>%
  ggplot(aes(x = mut, y = cluster, fill = er)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER', direction = 'horizontal')) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(vjust = 0.85))

p_surface_acc <- filter(dms, cluster %in% c('C1', 'C2')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density', show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_disulphide <- filter(dms, cluster %in% c('C1', 'C2'), !is.na(disulfide)) %>%
  count(cluster) %>%
  mutate(cluster = factor(cluster, levels = c('C1', 'C2')),
         cluster_int = as.integer(cluster)) %>%
  ggplot(aes(xmin = 0, xmax = n, ymin = as.integer(cluster) - 0.4, ymax = as.integer(cluster) + 0.4 , fill = cluster)) +
  geom_rect(show.legend = FALSE) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_disulphide.png'), interpolate = TRUE, xmin=-27, xmax=-7, ymin=1, ymax=2) +
  coord_fixed(ratio = 20, clip = 'off') +
  scale_y_continuous(breaks = c(1, 2), labels = c("<span style='color:#E41A1C'>C1</span>", "<span style='color:#377EB8'>C2</span>")) +
  scale_x_continuous(name = 'Disulphide Bonds', expand = expansion(0.01)) +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  theme(axis.text.y = element_markdown(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.22), 'npc'),
        plot.tag.position = c(0, 1))

p_ligand_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1.9), y = c(0, 0.9)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/gal4_cys_zinc.png'), interpolate = TRUE, xmin=0, xmax=0.9, ymin=0, ymax=0.9) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_cys_haem.png'), interpolate = TRUE, xmin=1, xmax=1.9, ymin=0, ymax=0.9) +
  labs(title = 'C1 Ligand Binding') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_aromatic_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1.9), y = c(0, 0.9)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/np_cys_aromatic.png'), interpolate = TRUE, xmin=0, xmax=0.9, ymin=0, ymax=0.9) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_aromatic.png'), interpolate = TRUE, xmin=1, xmax=1.9, ymin=0, ymax=0.9) +
  labs(title = 'C1/C2 Aromatic Interaction') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))

p1 <- p_profiles + size
p2 <- p_disulphide + size
p3 <- p_surface_acc + size
p4 <- p_ligand_ex + size
p5 <- p_aromatic_ex + size
 
figure4 <- multi_panel_figure(width = 200, height = 167, columns = 2, rows = 3,
                              panel_label_type = 'upper-alpha', row_spacing = 5, column_spacing = 5) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 3, column = 1) %>%
  fill_panel(p5, row = 3, column = 2)

ggsave('figures/4_figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
