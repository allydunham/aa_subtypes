#!/usr/bin/env Rscript
# Perform clustering on all positions at once
source('src/config.R')
source('src/subtypes.R')

dms_wide <- read_tsv('data/combined_mutational_scans.tsv')

### Make Dynamic Hclust Clusters ###
clusters <- make_dynamic_hclust_clusters(dms_wide, PC2:PC10, dist_method = 'manhattan', treecut_args = list(deepSplit=2))
dms_wide <- mutate(clusters$tbl, cluster = str_c('X', cluster))

### Analyse clusters ###
plots <- make_cluster_plots(dms_wide, cols = A:Y, chem_env_cols = within_10_0_A:within_10_0_Y, clusters = list(X=clusters))

cluster_occupancy <- group_by(dms_wide, cluster, wt) %>%
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
write_tsv(select(dms_wide, cluster, study, gene, position, wt), 'data/clustering/a.tsv')
root <- 'figures/2_clustering/hclust_profile_dynamic_all_positions'
dir.create(root)
save_plotlist(plots, root)
