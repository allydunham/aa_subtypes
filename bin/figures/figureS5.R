#!/usr/bin/env Rscript
# Produce figure S5 (Cluster consistency)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

cluster_order <- sort_clusters(unique(dms$cluster))

dupes <- mutate(dms, cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(desc(cluster)) %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(gene, position, wt) %>%
  filter(n() > 1) %>%
  summarise(clusters = str_c(cluster, collapse = ',')) %>%
  ungroup() %>%
  separate(clusters, c('cluster1', 'cluster2'), sep = ',') %>%
  mutate(same = if_else(cluster1 == cluster2, 'Match', 'Mismatch'))

p_overview <- ggplot(dupes, aes(x = same, fill=same)) +
  geom_bar(width=0.5) +
  scale_fill_manual(values = c(Match='cornflowerblue', Mismatch='firebrick2')) +
  guides(fill=FALSE) +
  labs(x = '', y = 'Count') +
  coord_flip() +
  theme(panel.grid.major.y = element_blank())

dupe_counts <- group_by(dupes, cluster1, cluster2) %>%
  tally() %>%
  ungroup() %>%
  mutate(cluster1 = factor(cluster1, levels = cluster_order),
         cluster2 = factor(cluster2, levels = cluster_order)) %>%
  complete(cluster1, cluster2) %>%
  mutate(n = ifelse(is.na(n), 0, n),
         wt = str_sub(cluster1, end = 1),
         match = if_else(cluster1 == cluster2, 'Match', 'Mismatch'),
         num1 = factor(str_sub(cluster1, 2), levels = c(1:10, 'P', 'O')),
         num2 = factor(str_sub(cluster2, 2), levels = c(1:10, 'P', 'O'))) %>%
  filter(str_sub(cluster2, end = 1) == wt,
         n > 0 | match == 'Match')

p_detail <- ggplot(dupe_counts, aes(x=num1, y=num2, size=n, colour=match)) +
  facet_wrap(~wt, scales = 'free') +
  geom_point() +
  coord_cartesian(clip = 'off') +
  scale_colour_manual(values = c(Match='cornflowerblue', Mismatch='firebrick2')) +
  scale_size_area() +
  scale_x_discrete() +
  scale_y_discrete() +
  labs(x='', y='') +
  theme(panel.grid.major.y = element_blank(),
        legend.position = 'top') +
  guides(size = guide_legend(title = ''), colour = guide_legend(title = ''))

figure <- multi_panel_figure(width = 183, height = c(160, 23), unit = 'mm', columns = 1) %>%
  fill_panel(p_detail, row = 1, column = 1) %>%
  fill_panel(p_overview, row = 2, column = 1)
ggsave('figures/4_figures/figureS5.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS5.png', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
