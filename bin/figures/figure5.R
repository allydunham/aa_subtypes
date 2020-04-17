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

charge_groups <- c(D1='Negative', D2='Polar', D3='Not Proline', E1='Negative', E2='Not Proline', E3='Polar',
                   K1='Positive', K2='Not Negative', K3='Not Proline', R1='Positive', R2='Not Proline', R3='Not Negative')

charge_profiles <- filter(full_characterisation$profiles, cluster %in% names(charge_groups)) %>%
  mutate(group = charge_groups[cluster]) %>%
  group_by(group, mut) %>%
  summarise(er = mean(er)) %>%
  ungroup() %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         group = as.factor(group))

charge_group_labs <- levels(charge_profiles$group)
charge_subtype_labs <- c(Negative="<span style = 'color:#e41a1c'>D1</span>, <span style = 'color:#e41a1c'>E1</span>",
                  Positive="<span style = 'color:#e41a1c'>K1</span>, <span style = 'color:#e41a1c'>R1</span>",
                  Polar="<span style = 'color:#377eb8'>D2</span>, <span style = 'color:#4daf4a'>E3</span>",
                  `Not Proline`="<span style = 'color:#4daf4a'>D3</span>, <span style = 'color:#377eb8'>E2</span>, <span style = 'color:#4daf4a'>K3</span>, <span style = 'color:#377eb8'>R2</span>",
                  `Not Negative`="<span style = 'color:#377eb8'>K2</span>, <span style = 'color:#4daf4a'>R3</span>")[charge_group_labs]

p_charge_profiles <- ggplot(charge_profiles, aes(x = mut, y = as.integer(group), fill = er)) +
  geom_raster(show.legend = FALSE) +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(charge_group_labs), labels = charge_group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(charge_subtype_labs), labels = charge_subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER')) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_charge_foldx <- filter(dms, cluster %in% names(charge_groups)) %>%
  mutate(group = charge_groups[cluster]) %>%
  ggplot(aes(x = group, y = electrostatics)) +
  geom_boxplot(show.legend = FALSE, fill = '#377eb8', outlier.shape = 20) +
  geom_hline(yintercept = 0, linetype = 'dotted', colour = 'black') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = '', y = expression('Mean Substitution Electrostatic'~Delta*Delta*'G (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))
er_legend <- plot_profiles('A1', legend = TRUE) %>% get_legend() %>% as_ggplot() + size

figure5 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size + 
  labs(tag = 'B') +
  annotation_custom(ggplotGrob(p_charge_profiles + size), xmin = 0, xmax = 1, ymin = 0.35, ymax = 0.68) +
  annotate('text', x = 0.5, y = 0.975, label = str_wrap('Charged subtype positions have the strongest electrostatic interactions', 80)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_asp_ionic.png'), interpolate = TRUE,
                    xmin = 0.025, xmax = 0.225, ymin = 0.725, ymax = 0.925) +
  annotate('text', x = 0.125, y = 0.7, label = str_wrap('Ionic interactions', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_ligand.png'), interpolate = TRUE,
                    xmin = 0.275, xmax = 0.475, ymin = 0.725, ymax = 0.925) +
  annotate('text', x = 0.375, y = 0.7, label = str_wrap('Ligand binding', 30)) + 
  annotation_custom(ggplotGrob(p_charge_foldx + size), xmin = 0.5, xmax = 1, ymin = 0.65, ymax = 0.925) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_sa.png'), interpolate = TRUE,
                    xmin = 0, xmax = 0.25, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.125, y = 0.3, label = str_wrap('Many polar positions are surface accessible', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_arg_not_proline.png'), interpolate = TRUE,
                    xmin = 0.375, xmax = 0.625, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.5, y = 0.3, label = str_wrap('Not proline positions are often in secondary structures or tight turns', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_arg_not_neg.png'), interpolate = TRUE,
                    xmin = 0.75, xmax = 1, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.875, y = 0.3, label = str_wrap('Positions near negative residues or ligands (e.g. DNA backbone) need to not be negative them selves', 30))


ggsave('figures/4_figures/figure5.pdf', figure5, width = 200, height = 200, units = 'mm')
ggsave('figures/4_figures/figure5.png', figure5, width = 200, height = 200, units = 'mm')
