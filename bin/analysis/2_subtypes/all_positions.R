#!/usr/bin/env Rscript
# Perform clustering on all positions at once
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

### Make Dynamic Hclust Clusters ###
clusters <- make_dynamic_hclust_clusters(dms, PC2:PC20, distance_method = 'cosine', treecut_args = list(deepSplit=0))
dms <- mutate(clusters$tbl,
              cluster = str_c('X', cluster) %>% relabel_outlier_clusters(),
              cluster = factor(cluster, levels = sort_clusters(cluster))) %>% # -1 as includes 0
  select(cluster, everything())

### Analyse clusters ###
n_clusters <- n_distinct(dms$cluster)
plots <- list()

plots$umap <- (ggplot(dms, aes(x=umap1, y=umap2, colour=wt)) +
                 facet_wrap(~cluster) +
                 scale_colour_manual(values = AA_COLOURS) +
                 geom_point() +
                 labs(x = 'UMAP1', y = 'UMAP2')) %>%
  labeled_plot(units='cm', height=25, width=25)

plots$tsne <- (ggplot(dms, aes(x=tSNE1, y=tSNE2, colour=wt)) +
                 facet_wrap(~cluster) +
                 geom_point() +
                 scale_colour_manual(values = AA_COLOURS)) %>%
  labeled_plot(units='cm', height=25, width=25)

plots$silhouette <- labeled_plot(plot_silhouette(dms, A:Y, 'cosine'),
                                 units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)


cluster_occupancy <- group_by(dms, cluster, wt) %>%
  tally() %>%
  mutate(rel_n = n / sum(n))
plots$cluster_occupancy <- filter(cluster_occupancy, !str_ends(cluster, '0')) %>%
  ggplot(aes(x = wt, y = cluster, fill = rel_n)) +
  geom_raster() +
  coord_equal() + 
  guides(fill = guide_colourbar(title = 'Proportion')) +
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = 1) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

### Save results ###
write_tsv(select(dms, cluster, study, gene, position, wt), 'data/subtypes/all_positions.tsv')
saveRDS(clusters, 'data/subtypes/all_positions.rds')

root <- 'figures/2_subtypes/all_positions'
dir.create(root)
save_plotlist(plots, root)
