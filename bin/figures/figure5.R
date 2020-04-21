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
  geom_tile(colour = 'grey') +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(group_labs), labels = group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(subtype_labs), labels = subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER', direction = 'horizontal')) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(vjust = 0.85))

p_foldx <- filter(dms, cluster %in% names(charge_groups)) %>%
  mutate(group = factor(charge_groups[cluster], levels = charge_group_order)) %>%
  ggplot(aes(x = group, y = electrostatics)) +
  geom_boxplot(show.legend = FALSE, fill = '#377eb8', outlier.shape = 20) +
  geom_hline(yintercept = 0, linetype = 'dotted', colour = 'black') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = '', y = expression('Mean Substitution Electrostatic'~Delta*Delta*'G (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

p_ionic_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_asp_ionic.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Ionic') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_negative_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_ligand.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Negative') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_positive_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/gal4_lys_dna.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Positive') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_polar_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_sa.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Polar') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_not_proline_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_arg_not_proline.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Not Proline') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_not_negative_ex <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_arg_not_neg.png'), interpolate = TRUE, xmin=0, xmax=1, ymin=0, ymax=1) +
  labs(title = 'Not Negative') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))

p1 <- p_profiles + size
p2 <- p_foldx + size
p3 <- p_ionic_ex + size
p4 <- p_polar_ex + size
p5 <- p_positive_ex + size
p6 <- p_negative_ex + size
p7 <- p_not_negative_ex + size
p8 <- p_not_proline_ex + size
 
figure5 <- multi_panel_figure(width = 200, height = 200, columns = 4, rows = 4,
                              panel_label_type = 'upper-alpha', row_spacing = 5, column_spacing = 5) %>%
  fill_panel(p1, row = 1:2, column = 1:4) %>%
  fill_panel(p2, row = 3, column = 1:2) %>%
  fill_panel(p3, row = 3, column = 3) %>%
  fill_panel(p4, row = 3, column = 4) %>%
  fill_panel(p5, row = 4, column = 1) %>%
  fill_panel(p6, row = 4, column = 2) %>%
  fill_panel(p7, row = 4, column = 3) %>%
  fill_panel(p8, row = 4, column = 4)

ggsave('figures/4_figures/figure5.pdf', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')
ggsave('figures/4_figures/figure5.png', figure5, width = figure_width(figure5), height = figure_height(figure5), units = 'mm')
