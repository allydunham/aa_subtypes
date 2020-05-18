#!/usr/bin/env Rscript
# Make and characterise clusters based on SIFT scores
source('src/config.R')
source('src/subtype_clustering.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')
sift <- read_tsv('data/long_combined_mutational_scans.tsv') %>%
  select(study, gene, position, wt, mut, sift) %>%
  filter(!mut == '*') %>%
  pivot_wider(names_from = mut, values_from = sift)

pca <- tibble_pca(sift, A:Y)
tsne <- tibble_tsne(sift, A:Y)
umap <- tibble_to_matrix(sift, A:Y) %>% umap(metric = 'manhattan')

sift <- bind_cols(tsne$tbl, as_tibble(pca$x), set_colnames(umap, c('umap1', 'umap2')) %>% as_tibble()) %>%
  left_join(.,
            select(dms, -matches('^.$'), -starts_with('PC'), -tSNE1, -tSNE2, -umap1, -umap2),
            by = c('study', 'gene', 'position', 'wt'))

# Filter permissive positions using SIFT > 0.05
permissive_positions <- tibble_to_matrix(sift, A:Y) %>%
  abs() %>%
  is_greater_than(0.05) %>%
  apply(1, all)

sift_permissive <- filter(sift, permissive_positions) %>%
  mutate(cluster = str_c(wt, CLUSTER_PERMISSIVE_CHAR))

sift <- filter(sift, !permissive_positions)

cluster_func <- function(tbl, ...){make_dynamic_hclust_clusters(tbl, PC2:PC20, distance_method = 'cosine', treecut_args = list('deepSplit'=2))}

clusters <- group_by(sift, wt) %>%
  group_map(cluster_func, keep = TRUE) %>%
  set_names(sapply(., function(x){first(x$tbl$wt)}))

sift <- map_dfr(clusters, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, cluster) %>% relabel_outlier_clusters()) %>%
  arrange(study, position) %>%
  bind_rows(sift_permissive)

### Diagnostic Plots ###
diagnostic_plots <- plot_cluster_diagnostics(sift, clusters, cols = A:Y)

### Calculate Profiles for all clusters ###
full_characterisation <- full_cluster_characterisation(sift)

# Make profiles with permissive/outlier clusters excluded
outlier_clusters <- filter(full_characterisation$profiles, str_detect(cluster, CLUSTER_PERMISSIVE_RE) | str_detect(cluster, CLUSTER_OUTLIER_RE)) %>%
  pull(cluster) %>%
  unique()
selective_characterisation <- full_cluster_characterisation(filter(sift, !cluster %in% outlier_clusters))
n_clusters_selective <- nrow(full_characterisation$summary)

### Plot all cluster characterisation ###
plots <- plot_cluster_characterisation(full_characterisation, selective_characterisation, clusters)

### Write Output ###
write_tsv(select(sift, cluster, study, gene, position, wt), 'data/subtypes/sift_scores.tsv')
saveRDS(clusters, file = 'data/subtypes/sift_scores.rds')
save_plotlist(diagnostic_plots, root = 'figures/2_subtypes/sift_scores', overwrite = 'all')
save_plotlist(plots, 'figures/2_subtypes/sift_scores', verbose = 2)

