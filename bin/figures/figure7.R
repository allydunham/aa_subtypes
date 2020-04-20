#!/usr/bin/env Rscript
# Produce figure 7 (Tyrosine Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

er_limits <- c(min(full_characterisation$profiles$er), -min(full_characterisation$profiles$er))

p_profiles <- filter(full_characterisation$profiles, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         cluster = add_markdown(cluster, cluster_number_colourmap(cluster))) %>%
  ggplot(aes(x = mut, y = cluster, fill = er)) +
  geom_raster(show.legend = TRUE) +
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

p_surface_acc <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density', show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

tyr_dssp <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  select(cluster, study, position, starts_with('ss_')) %>%
  pivot_longer(starts_with('ss_'), names_to = 'term', names_prefix = 'ss_', values_to = 'p') %>%
  group_by(cluster, study, position) %>%
  summarise(ss = DSSP_CLASSES_PLOTMATH[term[which.max(p)]]) %>%
  ungroup() %>%
  count(cluster, ss) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

p_dssp <- ggplot(tyr_dssp, aes(x = cluster, y = prop, fill = ss)) +
  geom_col() +
  scale_fill_brewer(labels = function(x){parse(text=x)}, type = 'qual', palette = 'Dark2') +
  coord_flip() +
  labs(y = 'Proportion', x = '') +
  guides(fill = guide_legend(title = '', reverse = TRUE, nrow = 1)) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')

p_ligand_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_tyr_haem.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Y1 Ligand Binding') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_pi_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/tp53_tyr_pi.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Y1 Pi Interaction') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_buried_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_tyr_buried.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Buried Y2') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_not_proline_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_tyr_not_proline.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Y4 Backbone Alignment') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))

p1 <- p_profiles + size
p2 <- p_surface_acc + size
p3 <- p_dssp + size
p4 <- p_ligand_ex + size
p5 <- p_pi_ex + size
p6 <- p_buried_ex + size
p7 <- p_not_proline_ex + size

figure7 <- multi_panel_figure(width = 200, height = 220, columns = 4, rows = 3,
                              panel_label_type = 'upper-alpha', row_spacing = 5, column_spacing = 5) %>%
  fill_panel(p1, row = 1, column = 1:4) %>%
  fill_panel(p2, row = 2, column = 1:2) %>%
  fill_panel(p3, row = 2, column = 3:4) %>%
  fill_panel(p4, row = 3, column = 1) %>%
  fill_panel(p5, row = 3, column = 2) %>%
  fill_panel(p6, row = 3, column = 3) %>%
  fill_panel(p7, row = 3, column = 4)
  
ggsave('figures/4_figures/figure7.pdf', figure7, width = figure_width(figure7), height = figure_height(figure7), units = 'mm')
ggsave('figures/4_figures/figure7.png', figure7, width = figure_width(figure7), height = figure_height(figure7), units = 'mm')
