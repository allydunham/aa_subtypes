#!/usr/bin/env Rscript
# Test methods for assigning new positions to clusters
source('src/config.R')

# Import Data
dms <- left_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by=c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, gene, position) %>%
  select(cluster:Y, PC1:PC20)

# Calculate "recluster" results, if these positions were assigned based on similarity to centroids
cluster_centers <- filter(dms, !str_detect(cluster, CLUSTER_OUTLIER_RE), !str_detect(cluster, CLUSTER_PERMISSIVE_RE)) %>%
  group_by(cluster) %>%
  summarise(across(.cols = c(A:Y, PC1:PC20), .fns = mean))

er_cluster_centers <- tibble_to_matrix(cluster_centers, A:Y, row_names = 'cluster')
pc_cluster_centers <- tibble_to_matrix(cluster_centers, PC2:PC20, row_names = 'cluster')

calc_perm_cluster <- function(wt, er){
  if (all(abs(er) < 0.4)){
    return(str_c(wt, 'P'))
  }
  return(NA)
}

calc_er_cluster <- function(wt, er){
  clus <- er_cluster_centers[str_starts(rownames(er_cluster_centers), wt), , drop = FALSE]
  m <- matrix(er, nrow = nrow(clus), ncol = 20, byrow = TRUE)
  d <- sqrt(rowSums((m - clus)^2))
  i <- which.min(d)
  return(str_c(d[i], '_', names(d)[i]))
}

calc_outlier_cutoff <- function(wt, pca){
  clus <- pc_cluster_centers[str_starts(rownames(pc_cluster_centers), wt), , drop = FALSE]
  m <- matrix(pca, nrow = nrow(clus), ncol = 19, byrow = TRUE)
  d <- acos(row_cosine_similarity(m, clus)) / pi
  if (all(d > 0.45)){
    return(str_c(wt, 'O'))
  }
  return(NA)
}

calc_outlier_margin <- function(wt, pca){
  clus <- pc_cluster_centers[str_starts(rownames(pc_cluster_centers), wt), , drop = FALSE]
  m <- matrix(pca, nrow = nrow(clus), ncol = 19, byrow = TRUE)
  d <- acos(row_cosine_similarity(m, clus)) / pi
  i <- which.min(d)
  diff <- d[-i] - d[i]
  if (length(diff) > 0) {
    if (min(diff) < 0.03) {
      return(str_c(wt, 'O'))
    }
  }
  return(NA)
}

calc_pca_cluster <- function(wt, pca){
  clus <- pc_cluster_centers[str_starts(rownames(pc_cluster_centers), wt), , drop = FALSE]
  m <- matrix(pca, nrow = nrow(clus), ncol = 19, byrow = TRUE)
  d <- acos(row_cosine_similarity(m, clus)) / pi
  i <- which.min(d)
  diff <- d[-i] - d[i]
  return(str_c(d[i], '_', ifelse(length(diff) > 0, min(diff), ''), '_', names(d)[i]))
}

combine_clusters <- function(base, perm=NA, cutoff=NA, margin=NA) {
  base[!is.na(cutoff)] <- cutoff[!is.na(cutoff)]
  base[!is.na(margin)] <- margin[!is.na(margin)]
  base[!is.na(perm)] <- perm[!is.na(perm)]
  return(base)
}

dms_recluster <- rowwise(dms) %>%
  mutate(perm_cluster = calc_perm_cluster(wt = wt, er = c(A, C, D, E, `F`, G, H, I, K, L, M, N, P, Q, R, S, `T`, V, W, Y)),
         er_cluster = calc_er_cluster(wt = wt, er = c(A, C, D, E, `F`, G, H, I, K, L, M, N, P, Q, R, S, `T`, V, W, Y)),
         pca_cluster = calc_pca_cluster(wt = wt, pca = c(PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14,
                                                         PC15, PC16, PC17, PC18, PC19, PC20)),
         outlier_cutoff = calc_outlier_cutoff(wt = wt, pca = c(PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14,
                                                                PC15, PC16, PC17, PC18, PC19, PC20)),
         outlier_margin = calc_outlier_margin(wt = wt, pca = c(PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14,
                                                               PC15, PC16, PC17, PC18, PC19, PC20))) %>%
  separate(er_cluster, c('er_dist', 'er_cluster'), sep = '_', convert = TRUE) %>%
  separate(pca_cluster, c('pca_dist', 'pca_margin', 'pca_cluster'), sep = '_', convert = TRUE) %>%
  mutate(er_perm_cluster = combine_clusters(er_cluster, perm=perm_cluster),
         pca_perm_cluster = combine_clusters(pca_cluster, perm=perm_cluster),
         pca_cutoff_cluster = combine_clusters(pca_cluster, cutoff=outlier_cutoff),
         pca_margin_cluster = combine_clusters(pca_cluster, margin=outlier_margin),
         pca_all_cluster = combine_clusters(pca_cluster, perm=perm_cluster, cutoff=outlier_cutoff, margin=outlier_margin)) %>%
  select(study, gene, position, wt, cluster, er_cluster, er_perm_cluster, pca_cluster, pca_margin_cluster,
         pca_cutoff_cluster, pca_perm_cluster, pca_all_cluster, er_dist, pca_dist, pca_margin)

# Analyse which method performs best, what distance suggests bad assignment, which clusters assign best/worst
recluster_summary <- select(dms_recluster, -er_dist, -pca_dist, -pca_margin) %>%
  rename(true = cluster) %>%
  pivot_longer(ends_with('cluster'), names_to = 'method', values_to = 'pred') %>%
  mutate(method = str_remove(method, '_cluster'))

calc_confusion <- function(tbl, ...){
  con <- sapply(unique(tbl$true), function(x){
    c(TP = sum(tbl$true == x & tbl$pred == x),
      TN = sum(tbl$true != x & tbl$pred != x),
      FP = sum(tbl$true != x & tbl$pred == x),
      FN = sum(tbl$true == x & tbl$pred != x)
  )})
  
  s <- rowSums(con)
  
  tibble(p_micro = s['TP'] / (s['TP'] + s['FP']),
         r_micro = s['TP'] / (s['TP'] + s['FN']),
         p_macro = sum(con['TP',] / (con['TP',] + con['FP',]), na.rm = TRUE) / ncol(con),
         r_macro = sum(con['TP',] / (con['TP',] + con['FN',]), na.rm = TRUE) / ncol(con))
}

metrics <- group_by(recluster_summary, method) %>%
  group_modify(calc_confusion) %>%
  ungroup() %>%
  mutate(f_micro = 2 * p_micro * r_micro / (p_micro + r_micro),
         f_macro = 2 * p_macro * r_macro / (p_macro + r_macro)) %>%
  left_join(group_by(recluster_summary, method) %>% summarise(accuracy = sum(true == pred) / n()), by = 'method') %>%
  left_join(filter(recluster_summary, !str_detect(pred, CLUSTER_OUTLIER_RE)) %>% group_by(method) %>% 
              summarise(accuracy_no_outlier = sum(true == pred) / n()), by = 'method') %>%
  arrange(accuracy, f_micro)

# Make plots
plots <- list()
plots$cutoff <- bind_rows(
  `ER` = select(dms_recluster, cluster, pred = er_cluster, dist = er_dist),
  `PCA` = select(dms_recluster, cluster, pred = pca_cluster, dist = pca_dist),
  `PCA Margin` = select(dms_recluster, cluster, pred = pca_cluster, dist = pca_margin),
  .id = 'metric') %>%
  mutate(true = cluster == pred) %>%
  ggplot(aes(x = dist, y = ..density.., colour = true)) +
  stat_density(geom = 'line', position = 'identity') +
  geom_vline(data = tibble(metric=c('PCA', 'PCA Margin'), x=c(0.45, 0.03)), mapping = aes(xintercept = x), linetype = 'dashed') +
  facet_wrap(~metric, scales = 'free', ncol = 1) +
  scale_colour_brewer(name = 'Match', type = 'qual', palette = 'Set1') +
  labs(x = 'Distance', y = 'Density')

plots$metrics <- (pivot_longer(metrics, -method, names_to = 'metric', values_to = 'value') %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  facet_wrap(~metric, strip.position = 'bottom', scales = 'free_x',
             labeller = labeller(metric = c(p_micro = 'Micro Precision', r_micro = 'Micro Recall',
                                            p_macro = 'Macro Precision', r_macro = 'Macro Recall',
                                            f_micro = 'Micro F1', f_macro = 'Macro F1',
                                            accuracy = 'Accuracy', accuracy_no_outlier = 'Accuracy (Excluding Outliers)'))) +
  geom_col(show.legend = FALSE, width = 0.7) +
  coord_flip() +
  lims(y = c(0,1)) +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  scale_x_discrete(labels = c(er = 'ER', er_perm = 'ER + Permissive', pca = 'PCA',
                              pca_outlier = 'PCA + Outliers', pca_perm = 'PCA + Permissive',
                              pca_outlier_perm = 'PCA + Outliers + Permissive')) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = 'outside',
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'lightgrey', linetype = 'dotted'))) %>%
  labeled_plot(units = 'cm', height = 20, width = 30)

save_plotlist(plots, root = 'figures/2_subtypes/cluster_assignment', overwrite = 'all')


