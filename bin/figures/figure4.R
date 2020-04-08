#!/usr/bin/env Rscript
# Produce figure 4 (Subtype Examples)
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
    mutate(cluster = add_markdown(cluster, cluster_colourmap(cluster)),
           mut = add_markdown(mut, AA_COLOURS)) %>%
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

### Example 1 - C examples ###
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
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

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

### Example 2 - D/E examples ###
p_charge_profiles <- plot_profiles(c('D1', 'D2', 'D3', 'E1', 'E2', 'E3'))

p_charge_foldx <- filter(dms, cluster %in% c('D1', 'D2', 'D3', 'E1', 'E2', 'E3'), !is.na(electrostatics)) %>%
  mutate(cluster_num = str_sub(cluster, -1)) %>%
  ggplot(aes(x = cluster, y = electrostatics, fill = cluster_num)) +
  geom_boxplot(show.legend = FALSE) +
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

### Example 3 - G1/G3 examples ###
p_gly_profiles <- plot_profiles(c('G1', 'G2', 'G3'))

p_gly_surface_acc <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density') +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_gly_foldx <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  select(cluster, entropy_sidechain, phi, psi, solvation_hydrophobic, solvation_polar, van_der_waals, van_der_waals_clashes) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'ddg') %>%
  mutate(term_pretty = c(FOLDX_TERMS, phi='Phi', psi='Psi')[term]) %>%
  ggplot(aes(x = cluster, y = ddg, fill = cluster)) +
  facet_wrap(~term_pretty, scales = 'free_x', nrow = 2) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill = FALSE) +
  labs(x = '', y = expression('Mean FoldX Prediction (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank())

p_gly_rama <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  ggplot(aes(x = phi, y = psi, colour = cluster)) +
  facet_wrap(~cluster) +
  geom_point() +
  coord_equal()

### Example 4 - Tyrosine ###
p_tyr_profiles <- plot_profiles(c('Y1', 'Y2', 'Y3', 'Y4', 'YP'))

p_tyr_surface_acc <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density') +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

tyr_dssp <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  select(cluster, ss_h, ss_e, ss_c) %>%
  pivot_longer(-cluster, names_to = 'term', names_prefix = 'ss_', values_to = 'p') %>%
  mutate(term = DSSP_CLASSES_PLOTMATH[term])

p_tyr_dssp <- ggplot(tyr_dssp, aes(x = cluster, y = p, fill = term)) +
  geom_boxplot() +
  scale_fill_brewer(labels = parse(text = sort(unique(tyr_dssp$term))), type = 'qual', palette = 'Dark2') +
  labs(y = 'Probability', x = '') +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank())

p_tyr_haem <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/cbs_tyr_haem.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_tyr_not_proline <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/pab1_tyr_not_proline.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_tyr_pi <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/tp53_tyr_pi.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_tyr_buried <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/aph3ii_tyr_buried.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_tyr <- blank_plot('Y Examples')

### Example 5 - Large/Small hydrophobics ###
large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
hydro_groups <- structure(rep(c('Large Aliphatic', 'Small Aliphatic', 'Aromatic'), each=3), names = c(large_hydrophobics, small_hydrophobics, aromatics))

p_hydro_profiles <- plot_profiles(names(hydro_groups))

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
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

p_hydro <- blank_plot('Large/Small/Aromatic Example')

### Assemble figure ###
size <- theme(text = element_text(size = 8))
er_legend <- plot_profiles('A1', legend = TRUE) %>% get_legend() %>% as_ggplot() + size

p1_prof <- p_cys_profiles + size + labs(tag = 'A')
p1_disulphide <- p_cys_disulphide + size
p1_disulphide_ex <- p_cys_disulphide_ex + size
p1_sa <- p_cys_surface_acc + size
p1_zinc <- p_cys_zinc_ex + size
p1_aromatic_1 <- p_cys_1_aromatic_ex + size
p1_aromatic_2 <- p_cys_2_aromatic_ex + size

p2 <- p_charge + labs(tag = 'B') + size

p3 <- p_tyr + labs(tag = 'C') + size

p4 <- p_hydro + labs(tag = 'D') + size

figure4 <- multi_panel_figure(width = 200, height = 200, columns = 14, rows = 14,
                              panel_label_type = 'none', row_spacing = 0.1, column_spacing = 0.1) %>%
  # Legend
  fill_panel(er_legend, row = 6:8, column = 13:14) %>%
  
  # Panel 1
  fill_panel(p1_prof, row = 1, column = 1:3) %>%
  fill_panel(p1_disulphide, row = 1, column = 4:5) %>%
  fill_panel(p1_disulphide_ex, row = 1, column = 6) %>%
  fill_panel(p1_sa, row = 2:3, column = 1:6) %>%
  fill_panel(p1_zinc, row = 4:6, column = 1:2) %>%
  fill_panel(p1_aromatic_1, row = 4:6, column = 3:4) %>%
  fill_panel(p1_aromatic_2, row = 4:6, column = 5:6) %>%
  
  # Panel 2 
  fill_panel(p2, row = 1:6, column = 7:12) %>%
  
  # Panel 3
  fill_panel(p3, row = 7:12, column = 1:6) %>%
  
  # Panel 4  
  fill_panel(p4, row = 7:12, column = 7:12)

ggsave('figures/4_figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
