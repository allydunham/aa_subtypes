#!/usr/bin/env Rscript
# Produce figure 6 (Large/Small Aliphatic Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

er_limits <- c(min(full_characterisation$profiles$er), -min(full_characterisation$profiles$er))

large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
hydro_groups <- structure(rep(c('Large Aliphatic', 'Small Aliphatic', 'Aromatic'), each=3), names = c(large_hydrophobics, small_hydrophobics, aromatics))

hydro_profiles <- filter(full_characterisation$profiles, cluster %in% names(hydro_groups)) %>%
  mutate(group = hydro_groups[cluster]) %>%
  group_by(group, mut) %>%
  summarise(er = mean(er)) %>%
  ungroup() %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         group = factor(group, levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic')))

group_labs <- add_markdown(levels(hydro_profiles$group), c(`Large Aliphatic`='#D95F02', `Small Aliphatic`='#7570B3', Aromatic='#1B9E77')) %>%
  str_c('**', ., '**')
subtype_labs <- map(levels(hydro_profiles$group),
                    ~names(hydro_groups[hydro_groups == .]) %>% 
                      add_markdown(., cluster_number_colourmap(.)) %>%
                      str_c(collapse = ', ')) %>%
  set_names(group_labs) %>% unlist()

p_profiles <- ggplot(hydro_profiles, aes(x = mut, y = as.integer(group), fill = er)) +
  geom_tile(colour = 'grey', size = 0.1) +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(group_labs), labels = group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(subtype_labs), labels = subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER', direction = 'horizontal', barheight = unit(2, 'mm'))) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.text.y.left = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, 'mm'),
        legend.position = 'bottom',
        legend.title = element_text(vjust = 0.85),
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(-3, 0, 0, 0, 'mm'),
        legend.background = element_blank())

p_foldx_entropy_sidechain <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = factor(hydro_groups[cluster], levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic'))) %>%
  ggplot(aes(x = as.integer(group), y = entropy_sidechain, fill = group)) +
  annotate('line', x = c(0, 3.75), y = 0, colour = 'black', linetype = 'dashed', size = 0.25) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ras_aliphatic_entropy.png'), interpolate = TRUE, xmin=4.5, xmax=8.5, ymin=-1, ymax=1) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, 0.5), limits = c(-1.5, 1.5)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 7)) +
  coord_fixed(ratio = 2, clip = 'off') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2', direction = -1) +
  labs(y = expression(Delta*Delta*'G (kj mol'^-1*')'), title = parse(text = FOLDX_TERMS_PLOTMATH['entropy_sidechain'])) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.length = unit(0.5, "mm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'mm'))

p_foldx_van_der_waals_clashes <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = factor(hydro_groups[cluster], levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic'))) %>%
  ggplot(aes(x = as.integer(group), y = van_der_waals_clashes, fill = group)) +
  annotate('line', x = c(0, 3.75), y = 0, colour = 'black', linetype = 'dashed', size = 0.25) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/adrb2_ala_small_hydro.png'), interpolate = TRUE, xmin=4.5, xmax=8.5, ymin=7.667, ymax=33.333) +
  lims(y = c(0, 40)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 7)) +
  coord_fixed(ratio = 0.15, clip = 'off') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2', direction = -1) +
  labs(y = expression(Delta*Delta*'G (kj mol'^-1*')'), title = parse(text = FOLDX_TERMS_PLOTMATH['van_der_waals_clashes'])) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.length = unit(0.5, "mm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'mm'))

p_foldx_solvation_hydrophobic <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = factor(hydro_groups[cluster], levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic'))) %>%
  ggplot(aes(x = as.integer(group), y = solvation_hydrophobic, fill = group)) +
  annotate('line', x = c(0, 3.75), y = 0, colour = 'black', linetype = 'dashed', size = 0.25) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ras_met_buried.png'), interpolate = TRUE, xmin=4.5, xmax=8.5, ymin=-2.667, ymax=2.667) +
  lims(y = c(-4, 4)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 7)) +
  coord_fixed(ratio = 0.75, clip = 'off') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2', direction = -1) +
  labs(y = expression(Delta*Delta*'G (kj mol'^-1*')'), title = parse(text = FOLDX_TERMS_PLOTMATH['solvation_hydrophobic'])) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.length = unit(0.5, "mm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'mm'))

p_foldx_van_der_waals <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = factor(hydro_groups[cluster], levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic'))) %>%
  ggplot(aes(x = as.integer(group), y = van_der_waals, fill = group)) +
  annotate('line', x = c(0, 3.75), y = 0, colour = 'black', linetype = 'dashed', size = 0.25) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_phe_pi.png'), interpolate = TRUE, xmin=4.5, xmax=8.5, ymin=-2, ymax=2) +
  scale_y_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 3)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 7)) +
  coord_fixed(clip = 'off') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2', direction = -1) +
  labs(y = expression(Delta*Delta*'G (kj mol'^-1*')'), title = parse(text = FOLDX_TERMS_PLOTMATH['van_der_waals'])) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.length = unit(0.5, "mm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'mm'))

p_sa <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = factor(hydro_groups[cluster], levels = c('Small Aliphatic', 'Large Aliphatic', 'Aromatic'))) %>%
  ggplot(aes(x = as.integer(group), y = all_atom_abs, fill = group)) +
  annotate('line', x = c(0, 3.75), y = 0, colour = 'black', linetype = 'dashed', size = 0.25) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_trp_surface.png'), interpolate = TRUE, xmin=4, xmax=8, ymin=41.667, ymax=208.333) +
  scale_y_continuous(limits = c(0, 250)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 7)) +
  coord_fixed(ratio = 0.024, clip = 'off') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2', direction = -1) +
  labs(y = 'All Atom Abs.', title = 'Surface Accessibility') +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        axis.ticks.length = unit(0.5, "mm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'mm'))

### Assemble figure ###
size <- theme(text = element_text(size = 6))

p1 <- p_profiles + size + labs(tag = 'A')
p2 <- p_foldx_entropy_sidechain + size + labs(tag = 'B')
p3 <- p_foldx_van_der_waals_clashes + size + labs(tag = 'C')
p4 <- p_foldx_van_der_waals + size + labs(tag = 'D')
p5 <- p_sa + size + labs(tag = 'E')

figure6 <- multi_panel_figure(width = 89, height = 89, columns = 2, rows = 3,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 3, column = 1) %>%
  fill_panel(p5, row = 3, column = 2)
  
ggsave('figures/4_figures/figure6.pdf', figure6, width = figure_width(figure6), height = figure_height(figure6), units = 'mm')
ggsave('figures/4_figures/figure6.png', figure6, width = figure_width(figure6), height = figure_height(figure6), units = 'mm')
ggsave('figures/4_figures/figure6.tiff', figure6, width = figure_width(figure6), height = figure_height(figure6), units = 'mm')
ggsave('figures/4_figures/figure6.eps', figure6, width = figure_width(figure6), height = figure_height(figure6), units = 'mm', device=cairo_ps, fallback_resolution = 600)
