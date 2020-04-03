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
most_selective_subtypes <- group_by(full_characterisation$profiles, cluster) %>%
  summarise(mean_er = mean(er)) %>%
  mutate(aa = str_sub(cluster, end = 1)) %>%
  group_by(aa) %>%
  filter(!mean_er > min(mean_er))

cor_subtype_sets <- list(small_aliphatic = c('G1', 'G2', 'A5', 'G3', 'S1'),
                         not_proline = c('I3', 'Y4', 'T2', 'E2', 'D3', 'L6', 'N2', 'R2', 'K3', 'A3', 'V5', 'S2', 'Q2', 'M2'),
                         positive_polar = c('Q4', 'K5', 'R3', 'K2', 'R5', 'K1', 'R1'),
                         aromatic = c('Y3', 'Y1', 'H1', 'W1', 'F2'),
                         aliphatic = c('C2', 'A2', 'G7', 'S3', 'V1', 'L3', 'I1', 'P4', 'L4'),
                         large_aliphatic = c('P4', 'L4', 'T4', 'L5', 'L1', 'F1', 'Y2', 'V3', 'V2', 'M1', 'I2', 'L2'),
                         not_aromatic = c('T6', 'P2', 'T1', 'A1', 'P3', 'A4', 'T3', 'R4'),
                         negative_polar = c('T5', 'Q3', 'N1', 'Q1', 'E3', 'L7', 'G4', 'D2', 'D1', 'E1'))

### Panel 1 - Schematic ###
## Subparts
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

(pivot_wider(full_characterisation$profiles, cluster, names_from = 'mut', values_from = 'er') %>%
    filter(str_starts(cluster, 'A')) %>%
    plot_profile_dendogram(A:Y) +
    guides(colour = FALSE)) %>%
  ggsave('figures/4_figures/parts/figure3_cluster_schematic_dend.pdf', ., width = 10, height = 7, units = 'cm')

## Full Figure
p_schematic <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/parts/figure3_cluster_schematic.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

### Panel 2 - Cluster Sizes ###
subtype_size <- group_by(dms, cluster, wt) %>%
  summarise(size = n()) %>%
  group_by(wt) %>%
  mutate(n_aa = sum(size)) %>%
  ungroup() %>%
  mutate(freq = size / n_aa,
         cluster_num = str_sub(cluster, -1),
         wt = add_markdown(wt, AA_COLOURS),
         cluster_num = factor(cluster_num, levels = c('O', 'P', 8:1)))

aa_labs <- levels(subtype_size$wt)
n_aa_labs <- structure(subtype_size$n_aa, names = as.character(subtype_size$wt))[aa_labs]

p_sizes <- ggplot(subtype_size, aes(y = as.integer(wt), x = freq, fill = cluster_num)) +
  geom_col() +
  scale_y_continuous(breaks = 1:length(aa_labs), labels = aa_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(n_aa_labs), labels = n_aa_labs), name = ' Total Positions') +
  scale_fill_brewer(type = 'qual', palette = 'Paired') +
  guides(fill = guide_legend(title = 'Subtype', reverse = TRUE)) +
  labs(x = 'Subtype Frequency') +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_markdown(),
        panel.grid.major.y = element_blank(),
        legend.position = 'top')

### Panel 3 - Correlation plot/dendrogram with labels ###
## Subparts
# Correlation heatmap/dendrogram
cors <- filter(dms, !str_detect(cluster, CLUSTER_OUTLIER_RE), !str_detect(cluster, CLUSTER_PERMISSIVE_RE)) %>%
  cluster_profile_correlation(A:Y)

hc <- select(cors, cluster1, cluster2, cor) %>%
  pivot_wider(id_cols = cluster1, names_from = cluster2, values_from = cor) %>%
  tibble_to_matrix(-cluster1, row_names = 'cluster1') %>%
  subtract(1, .) %>%
  as.dist() %>%
  hclust()
dend_data <- dendro_data(hc)
branches <- dend_data$segments
leaves <- dend_data$labels %>%
  mutate(wt = str_sub(label, end = 1))

cors <- mutate(cors,
               cluster1 = factor(cluster1, levels = leaves$label),
               cluster2 = factor(cluster2, levels = leaves$label))

p_cor_heatmap <- ggplot(cors, aes(x=cluster1, y=cluster2, fill=cor)) +
  geom_tile() +
  scale_fill_distiller(type = ER_COR_COLOURS$type, palette = ER_COR_COLOURS$palette, direction = ER_COR_COLOURS$direction,
                       limits = c(-1, 1)) +
  coord_fixed() +
  guides(fill = guide_colourbar(title = 'Pearson\nCorrelation')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())
ggsave('figures/4_figures/parts/figure3_cor_heatmap.pdf', p_cor_heatmap, width = 12, height = 12, units = 'cm')

p_cor_dend <- ggplot() +
  geom_segment(data = branches, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_point(data = leaves, aes(x=x, y=y, colour = wt), shape=19) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  scale_colour_manual(values = AA_COLOURS) +
  guides(colour = FALSE) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())
ggsave('figures/4_figures/parts/figure3_cor_dendrogram.pdf', p_cor_dend, width = 12, height = 3, units = 'cm')

# Profiles of key groups
plot_profile_block <- function(set){
  filter(full_characterisation$profiles, cluster %in% set) %>%
    mutate(cluster = add_markdown(cluster, cluster_colourmap(cluster), order = levels(leaves$label)[levels(leaves$label) %in% set]),
           mut = add_markdown(mut, AA_COLOURS)) %>%
    ggplot(aes(x = mut, y = cluster, fill = er)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                         limits = c(-1, 1)) +
    guides(fill = guide_colourbar(title = 'Mean ER')) +
    theme(text = element_text(size = 8),
          axis.ticks = element_blank(),
          axis.text.y = element_markdown(),
          axis.text.x = element_markdown(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          panel.grid.major.y = element_blank())
}

ggsave('figures/4_figures/parts/figure3_cor_set_small_aliphatic.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$small_aliphatic))
ggsave('figures/4_figures/parts/figure3_cor_set_not_proline.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$not_proline))
ggsave('figures/4_figures/parts/figure3_cor_set_positive.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$positive_polar))
ggsave('figures/4_figures/parts/figure3_cor_set_aromatic.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$aromatic))
ggsave('figures/4_figures/parts/figure3_cor_set_aliphatic.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$aliphatic))
ggsave('figures/4_figures/parts/figure3_cor_set_larger_aliphatic.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$large_aliphatic))
ggsave('figures/4_figures/parts/figure3_cor_set_not_aromatic.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$not_aromatic))
ggsave('figures/4_figures/parts/figure3_cor_set_negative.pdf', width = 12, height = 6, units = 'cm',
       plot = plot_profile_block(cor_subtype_sets$negative_polar))

## Main plot
p_heatmap <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/parts/figure3_cor.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

### Panel 4 - Clusters mapped to UMAP ###
classify_cluster <- function(x){
  out <- rep(NA, length(x))
  out[str_detect(x, CLUSTER_PERMISSIVE_RE)] <- 'Permissive'
  out[x %in% cor_subtype_sets$not_proline] <- 'Not Proline'
  out[x %in% cor_subtype_sets$large_aliphatic] <- 'Large Aliphatic'
  out[x %in% cor_subtype_sets$small_aliphatic] <- 'Small Aliphatic'
  out[x %in% cor_subtype_sets$aliphatic] <- 'Aliphatic'
  out[x %in% cor_subtype_sets$positive_polar] <- 'Polar (Positive)'
  out[x %in% cor_subtype_sets$negative_polar] <- 'Polar (Negative)'
  out[x %in% cor_subtype_sets$aromatic] <- 'Aromatic'
  out[x %in% cor_subtype_sets$not_aromatic] <- 'Not Aromatic'
  return(out)
}

cor_set_colours <- c(`Small Aliphatic` = '#ff7f00', Aliphatic = '#ffff33', `Large Aliphatic` = '#a65628',
                     `Not Proline` = '#4daf4a', Permissive = '#999999',
                     `Polar (Positive)` = '#984ea3', `Polar (Negative)` = '#f781bf',
                     Aromatic = '#377eb8', `Not Aromatic` = '#e41a1c') # Manual assignment of colourbrewer2 Set1 colours

cluster_dms <- mutate(dms, cluster_type = classify_cluster(cluster)) %>%
  select(gene, position, wt, umap1, umap2, cluster, cluster_type) %>%
  drop_na(cluster_type)

p_umap <- ggplot() +
  geom_point(data = dms, mapping = aes(x = umap1, y = umap2), colour = 'grey90', shape = 20) +
  geom_point(data = cluster_dms, mapping = aes(x = umap1, y = umap2, colour = cluster_type)) +
  scale_color_manual(values = cor_set_colours) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_legend(title = '')) +
  theme(legend.position = 'top',
        legend.title = element_blank())

### Panel 5 - Subtype frequencies ###
freq_summary <- group_by(full_characterisation$summary, aa) %>%
  mutate(freq = n / sum(n)) %>%
  summarise(Permissive = freq[which(cluster == str_c(aa, 'P'))],
            `Not Proline` = max(freq[which(cluster %in% cor_subtype_sets$not_proline)], 0),
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

### Panel 6 - Profile of most selective subtype ###
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
p2 <- p_sizes + labs(tag = 'B') + size
p3 <- p_heatmap + labs(tag = 'C') + size
p4 <- p_umap + labs(tag = 'D') + size
p5 <- p_subtype_freqs + labs(tag = 'E') + size
p6 <- p_selective_profiles + labs(tag = 'F') + size

figure3 <- multi_panel_figure(width = 300, height = 300, columns = 9, rows = 3,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1:3) %>%
  fill_panel(p2, row = 2, column = 1:3) %>%
  fill_panel(p3, row = 1:2, column = 4:9) %>%
  fill_panel(p4, row = 3, column = 1:3) %>%
  fill_panel(p5, row = 3, column = 4:6) %>%
  fill_panel(p6, row = 3, column = 7:9)
ggsave('figures/4_figures/figure3.pdf', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm')
ggsave('figures/4_figures/figure3.png', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm')
