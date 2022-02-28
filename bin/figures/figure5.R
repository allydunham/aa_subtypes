#!/usr/bin/env Rscript
# Produce figure 5 (Charged Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

er_limits <- c(min(full_characterisation$profiles$er), -min(full_characterisation$profiles$er))

charge_groups <- c(D1='Negative', D2='Polar', D3='Not Proline', E1='Negative', E2='Not Proline', E3='Polar',
                   K1='Positive', K2='Not Negative', K3='Not Proline', R1='Positive', R2='Not Proline', R3='Not Negative')

charge_group_order <- mutate(dms, g = charge_groups[cluster]) %>%
  drop_na(g) %>%
  group_by(g) %>%
  summarise(f = mean(electrostatics, na.rm = TRUE)) %>%
  arrange(f) %>%
  pull(g)

charge_profiles <- filter(full_characterisation$profiles, cluster %in% names(charge_groups)) %>%
  mutate(group = charge_groups[cluster]) %>%
  group_by(group, mut) %>%
  summarise(er = mean(er)) %>%
  ungroup() %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         group = factor(group, levels = charge_group_order))

group_labs <- levels(charge_profiles$group)
subtype_labs <- c(Negative="<span style = 'color:#e41a1c'>D1</span>, <span style = 'color:#e41a1c'>E1</span>",
                  Positive="<span style = 'color:#e41a1c'>K1</span>, <span style = 'color:#e41a1c'>R1</span>",
                  Polar="<span style = 'color:#377eb8'>D2</span>, <span style = 'color:#4daf4a'>E3</span>",
                  `Not Proline`="<span style = 'color:#4daf4a'>D3</span>, <span style = 'color:#377eb8'>E2</span>, <span style = 'color:#4daf4a'>K3</span>, <span style = 'color:#377eb8'>R2</span>",
                  `Not Negative`="<span style = 'color:#377eb8'>K2</span>, <span style = 'color:#4daf4a'>R3</span>")[group_labs]

p_profiles <- ggplot(charge_profiles, aes(x = mut, y = as.integer(group), fill = er)) +
  geom_tile(colour = 'grey', size = 0.1) +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(group_labs), labels = group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(subtype_labs), labels = subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER', direction = 'horizontal', barheight = unit(2, 'mm'))) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, 'mm'),
        legend.position = 'bottom',
        legend.title = element_text(vjust = 0.85),
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(-3, 0, 0, 0, 'mm'),
        legend.background = element_blank())

p_foldx <- filter(dms, cluster %in% names(charge_groups)) %>%
  mutate(group = factor(charge_groups[cluster], levels = charge_group_order)) %>%
  ggplot(aes(x = group, y = electrostatics)) +
  geom_boxplot(show.legend = FALSE, fill = '#377eb8', outlier.shape = 20, outlier.size = 0.25, lwd = 0.1) +
  geom_hline(yintercept = 0, linetype = 'dotted', colour = 'black', size = 0.25) +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = '', y = expression('Electrostatic'~Delta*Delta*'G (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.5, "mm"),
        axis.title.x = element_text(hjust = 1))

p_ionic_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_asp_ionic.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Ionic') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

p_negative_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_ligand.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Negative') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

p_positive_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/gal4_lys_dna.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Positive') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

p_polar_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_sa.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Polar') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

p_not_proline_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_arg_not_proline.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Not Proline') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

p_not_negative_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_arg_not_neg.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Not Negative') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(margin = margin(1, 1, 1, 1, 'mm')),
        plot.margin = unit(c(2, 2, 2, 2), 'mm'))

### Assemble figure ###
size <- theme(text = element_text(size = 7))

p1 <- p_profiles + labs(tag = 'A') + size
p2 <- p_foldx + labs(tag = 'B') + size
p3 <- p_polar_ex + labs(tag = 'C') + size
p4 <- p_ionic_ex + labs(tag = 'D') + size
p5 <- p_positive_ex + labs(tag = 'E') + size
p6 <- p_negative_ex + labs(tag = 'F') + size
p7 <- p_not_negative_ex + labs(tag = 'G') + size
p8 <- p_not_proline_ex + labs(tag = 'H') + size
 
figure5 <- multi_panel_figure(width = 89, height = 89, columns = 4, rows = 3,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:4) %>%
  fill_panel(p2, row = 2, column = 1:2) %>%
  fill_panel(p3, row = 2, column = 3) %>%
  fill_panel(p4, row = 2, column = 4) %>%
  fill_panel(p5, row = 3, column = 1) %>%
  fill_panel(p6, row = 3, column = 2) %>%
  fill_panel(p7, row = 3, column = 3) %>%
  fill_panel(p8, row = 3, column = 4)

ggsave('figures/4_figures/figure5.pdf', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')
ggsave('figures/4_figures/figure5.png', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')
ggsave('figures/4_figures/figure5.tiff', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')
ggsave('figures/4_figures/figure5.eps', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm', device=cairo_ps, fallback_resolution = 600)

