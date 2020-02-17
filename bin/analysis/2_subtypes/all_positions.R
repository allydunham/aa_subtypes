#!/usr/bin/env Rscript
# Perform clustering on all positions at once
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

### Make Dynamic Hclust Clusters ###
clusters <- make_dynamic_hclust_clusters(dms, PC2:PC20, distance_method = 'cosine', treecut_args = list(deepSplit=0))
dms <- mutate(clusters$tbl, cluster = factor(str_c('X', cluster), levels = str_c('X', c(1:(n_distinct(cluster) - 1), 0)))) %>% # -1 as includes 0
  select(cluster, everything())

### Analyse clusters ###
n_clusters <- n_distinct(dms$cluster)

plots <- list()
plots$umap <- labeled_plot(plot_cluster_umap(dms), units = 'cm', height = 20, width = 20)
plots$tsne <- labeled_plot(plot_cluster_tsne(dms), units = 'cm', height = 20, width = 20)
plots$silhouette_global <- labeled_plot(plot_silhouette(dms, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
plots$silhouette_per_aa <- labeled_plot(plot_per_aa_silhouette(dms, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
plots$silhouette_per_aa_cosine <- labeled_plot(plot_per_aa_silhouette(dms, A:Y, 'cosine'), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)

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
write_tsv(select(dms, cluster, study, gene, position, wt), 'data/clustering/a.tsv')
root <- 'figures/2_clustering/hclust_profile_dynamic_all_positions'
dir.create(root)
save_plotlist(plots, root)
