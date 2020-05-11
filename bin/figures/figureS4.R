#!/usr/bin/env Rscript
# Produce figure S4 (COnsistency of positions in UMAP space)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')
dms <- left_join(dms, count(dms, gene, position), by = c('gene', 'position'))

repeated <- filter(dms, n > 1) %>%
  select(study, gene, position, umap1, umap2) %>%
  group_by(gene, position) %>%
  summarise(umap1_1 = umap1[1],
            umap1_2 = umap1[2],
            umap2_1 = umap2[1],
            umap2_2 = umap2[2])

p_umap <- ggplot() +
  geom_point(data = dms, mapping = aes(x = umap1, y = umap2), colour = 'grey90', shape = 20, size = 0.8) +
  geom_segment(data = repeated, aes(x = umap1_1, y = umap2_1, xend = umap1_2, yend = umap2_2)) +
  geom_point(data = filter(dms, n > 1), mapping = aes(x = umap1, y = umap2, colour = gene)) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  guides(colour = guide_legend(title = '')) +
  labs(x = 'UMAP1', y = 'UMAP2')

# Distribution of distances n > 1/n == 1
distances <- mutate(dms, gene_pos = str_c(gene, '_', position)) %>%
  tibble_to_matrix(umap1, umap2, row_names = 'gene_pos') %>%
  dist() %>%
  as.matrix()
distances[upper.tri(distances, diag = TRUE)] <- NA
distances <- as_tibble(distances, rownames = 'gene_pos1') %>%
  pivot_longer(-gene_pos1, names_to = 'gene_pos2', values_to = 'dist') %>%
  drop_na() %>%
  mutate(rep = ifelse(gene_pos1 == gene_pos2, 'Repeated Position', 'Background'))

# t.test
# wilcox.test(x = filter(distances, rep == 'Repeated Position')$dist, y = filter(distances, rep == 'Background')$dist, alternative = 'less')

p_dists <- ggplot(distances, aes(x = dist, y = ..scaled.., colour = rep)) +
  stat_density(geom = 'line', position = 'identity') +
  labs(x = 'UMAP Space Euclidean Distance', y = 'Scaled Density') + 
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  guides(colour = guide_legend(title = ''))
  
figure <- multi_panel_figure(width = 183, height = c(89, 89), unit = 'mm', columns = 1) %>%
  fill_panel(p_umap, row = 1, column = 1) %>%
  fill_panel(p_dists, row = 2, column = 1)
ggsave('figures/4_figures/figureS4.pdf', figure, width = 183, height = 185, units = 'mm')
ggsave('figures/4_figures/figureS4.png', figure, width = 183, height = 185, units = 'mm')

