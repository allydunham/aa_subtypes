#!/usr/bin/env Rscript
# Produce figure 3 (Subtype Overview)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

# Subtype Categories
not_proline_subtypes <- c('A3', 'D3', 'E2', 'G4', 'I3', 'K3', 'L6', 'M2', 'N2', 'Q2', 'R2', 'S2', 'T2', 'V5', 'Y4')
most_selective_subtypes <- group_by(full_characterisation$profiles, cluster) %>%
  summarise(mean_er = mean(er)) %>%
  mutate(aa = str_sub(cluster, end = 1)) %>%
  group_by(aa) %>%
  filter(!mean_er > min(mean_er))

large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
polar <- c('D1', 'D2', 'E1', 'K1', 'N1', 'Q1', 'Q3', 'R1', 'S1', 'T1')

### Panel 1 - Schematic? ###
position_profs <- filter(dms, wt %in% c('A', 'C', 'D', 'W', 'Y')) %>%
  select(study, cluster, position, wt, A:Y) %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
  mutate(er = clamp(er, 1, -1))

p_initial_profiles <- ggplot(position_profs, aes(x = str_c(study, position), y = mut, fill = er)) +
  geom_raster() +
  facet_wrap(~wt, nrow = 1, scales = 'free_x') +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                       limits = c(min(position_profs$er), -min(position_profs$er))) +
  guides(fill = FALSE) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 6))
ggsave('figures/4_figures/parts/figure3_cluster_schematic_initial_profiles.pdf', p_initial_profiles, width = 6, height = 2, units = 'cm')

p_a_permissive_profiles <- filter(position_profs, wt == 'A') %>%
  mutate(cat = ifelse(cluster == 'AP', 'AP', 'Rest')) %>%
  ggplot(aes(x = str_c(study, position), y = mut, fill = er)) +
  geom_raster() +
  facet_wrap(~cat, nrow = 1, scales = 'free_x') +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                       limits = c(min(position_profs$er), -min(position_profs$er))) +
  guides(fill = FALSE) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 6))
ggsave('figures/4_figures/parts/figure3_cluster_schematic_permissive_profs.pdf', p_a_permissive_profiles, width = 3, height = 2, units = 'cm')

p_schematic <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/parts/figure3_cluster_schematic.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

### Panel 2 - Correlation plot/dendrogram with labels ###
p_heatmap <- blank_plot('Correlation Heatmap')

### Panel 3 - Clusters mapped to UMAP ###
classify_cluster <- function(x){
  out <- rep(NA, length(x))
  out[str_detect(x, CLUSTER_PERMISSIVE_RE)] <- 'Permissive'
  out[x %in% not_proline_subtypes] <- 'Not Proline'
  out[x %in% large_hydrophobics] <- 'Large Hydrophobic'
  out[x %in% small_hydrophobics] <- 'Small Hydrophobic'
  out[x %in% polar] <- 'Polar'
  out[x %in% aromatics] <- 'Aromatic'
  return(out)
}

cluster_dms <- mutate(dms, cluster_type = classify_cluster(cluster)) %>%
  select(gene, position, wt, umap1, umap2, cluster, cluster_type) %>%
  drop_na(cluster_type)

p_umap <- ggplot() +
  geom_point(data = dms, mapping = aes(x = umap1, y = umap2), colour = 'grey90', shape = 20) +
  geom_point(data = cluster_dms, mapping = aes(x = umap1, y = umap2, colour = cluster_type)) +
  scale_color_brewer(type = 'qual', palette = 'Dark2') +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_legend(title = 'Subtype tolerates:'))

### Panel 4 - Subtype frequencies ###
freq_summary <- group_by(full_characterisation$summary, aa) %>%
  mutate(freq = n / sum(n)) %>%
  summarise(Permissive = freq[which(cluster == str_c(aa, 'P'))],
            `Not Proline` = max(freq[which(cluster %in% not_proline_subtypes)], 0),
            `Most Selective` = freq[which(cluster %in% most_selective_subtypes$cluster)],
            Other = 1 - (Permissive + `Not Proline` + `Most Selective`)) %>%
  pivot_longer(-aa, names_to = 'type', values_to = 'freq') %>%
  mutate(type = factor(type, levels = c('Permissive', 'Not Proline', 'Other', 'Most Selective'))) %>%
  left_join(select(most_selective_subtypes, aa, mean_er), by = 'aa') %>%
  mutate(aa = add_markdown(aa, AA_COLOURS))

p_subtype_freqs <- ggplot(freq_summary) +
  geom_col(aes(x = freq, y = aa, fill = type)) +
  scale_fill_brewer(type = 'qual', palette = 'Paired') + 
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = 'Frequency') +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(),
        axis.ticks.y = element_blank())

### Panel 5 - Profile of most selective subtype ###
most_selective_profiles <- filter(full_characterisation$profiles, cluster %in% most_selective_subtypes$cluster) %>%
  mutate(mut = add_markdown(mut, colour = AA_COLOURS),
         cluster = factor(cluster, levels = rev(sort(unique(cluster)))))

p_selective_profiles <- ggplot(most_selective_profiles, aes(y = cluster, x = mut, fill = er)) +
  geom_raster() +
  facet_wrap(~cluster, ncol = 1, scales = 'free_y', strip.position = 'left') +
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = -1, limits = c(-0.8, max(most_selective_profiles$er))) +
  labs(x = 'Substitution', y = '') +
  guides(fill = guide_colourbar(title = 'ER')) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_markdown(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.1, 'lines'),
        panel.grid.major.y = element_blank())

### Assemble figure ###
size <- theme(text = element_text(size = 8))
p1 <- p_schematic + labs(tag = 'A') + size
p2 <- p_heatmap + labs(tag = 'B') + size
p3 <- p_umap + labs(tag = 'C') + size
p4 <- p_subtype_freqs + labs(tag = 'D') + size
p5 <- p_selective_profiles + labs(tag = 'E') + size

figure3 <- multi_panel_figure(width = 300, height = 300, columns = 9, rows = 3,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1:2, column = 1:3) %>%
  fill_panel(p2, row = 1:2, column = 4:9) %>%
  fill_panel(p3, row = 3, column = 1:3) %>%
  fill_panel(p4, row = 3, column = 4:6) %>%
  fill_panel(p5, row = 3, column = 7:9)
ggsave('figures/4_figures/figure3.pdf', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm')
ggsave('figures/4_figures/figure3.png', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm')
