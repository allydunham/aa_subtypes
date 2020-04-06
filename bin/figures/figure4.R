#!/usr/bin/env Rscript
# Produce figure 4 (Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

### Panel 1 - C examples ###
p_cys_surface_acc <- filter(dms, cluster %in% c('C1', 'C2')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density') +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_cys_disulphide <- filter(dms, cluster %in% c('C1', 'C2'), !is.na(disulfide)) %>%
  ggplot(aes(x = cluster, fill = cluster)) +
  geom_bar() +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  labs(y = 'Count', x = '') +
  theme(legend.title = element_blank())

p_cys_disulphide_ex <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/ccr5_cys_disulphide.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_cys_2_aromatic_ex <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/ccr5_cys_aromatic.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_cys_1_aromatic_ex <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/np_cys_aromatic.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_cys_zinc_ex <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/gal4_cys_zinc.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_cys <- blank_plot('C examples')

### Panel 2 - D/E examples ###
p_charge_foldx <- filter(dms, cluster %in% c('D1', 'D2', 'D3', 'E1', 'E2', 'E3'), !is.na(electrostatics)) %>%
  mutate(cluster_num = str_sub(cluster, -1)) %>%
  ggplot(aes(x = wt, y = electrostatics, fill = cluster_num)) +
  geom_boxplot() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  coord_flip() +
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = '', y = expression('Mean FoldX Electrostatic'~Delta*'G (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

p_asp_ionic <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/cbs_asp_ionic.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_asp_ligand <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/tem1_asp_ligand.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_charge <- blank_plot('D/E (R/K) Examples')

### Panel 3 - G1/G3 examples ###
p_gly_foldx <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  ggplot(aes(x = phi, y = psi, colour = cluster)) +
  geom_point() +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  guides(colour = guide_legend(title = 'Subtype')) +
  labs(x = 'Phi', y = 'Psi') +
  theme()

p_gly <- blank_plot('G1/G3 Example')

### Panel 4 - Large/Small hydrophobics ###
large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
hydro_groups <- structure(rep(c('Large Aliphatic', 'Small Aliphatic', 'Aromatic'), each=3), names = c(large_hydrophobics, small_hydrophobics, aromatics))

p_hydro_foldx <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  select(cluster, entropy_mainchain, entropy_sidechain, solvation_polar, solvation_hydrophobic,
         van_der_waals, van_der_waals_clashes, torsional_clash, backbone_clash) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'ddg') %>%
  mutate(group = hydro_groups[cluster],
         term_pretty = FOLDX_TERMS[term]) %>%
  ggplot(aes(x = group, y = ddg, fill = group)) +
  facet_wrap(~term_pretty, scales = 'free_x', nrow = 2) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill = FALSE) +
  labs(x = '', y = expression('Mean FoldX Prediction (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank())

p_hydro <- blank_plot('Large/Small/Aromatic Example')

### Assemble figure ###
size <- theme(text = element_text(size = 8))
p1 <- p_cys + labs(tag = 'A') + size
p2 <- p_charge + labs(tag = 'B') + size
p3 <- p_gly + labs(tag = 'C') + size
p4 <- p_hydro + labs(tag = 'D') + size

figure4 <- multi_panel_figure(width = 200, height = 200, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)
ggsave('figures/4_figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
