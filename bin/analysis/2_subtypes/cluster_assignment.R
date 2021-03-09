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

calc_pca_cluster <- function(wt, pca){
  clus <- pc_cluster_centers[str_starts(rownames(pc_cluster_centers), wt), , drop = FALSE]
  m <- matrix(pca, nrow = nrow(clus), ncol = 19, byrow = TRUE)
  d <- acos(row_cosine_similarity(m, clus)) / pi
  i <- which.min(d)
  return(str_c(d[i], '_', names(d)[i]))
}

dms_recluster <- rowwise(dms) %>%
  mutate(perm_cluster = calc_perm_cluster(wt = wt, er = c(A, C, D, E, `F`, G, H, I, K, L, M, N, P, Q, R, S, `T`, V, W, Y)),
         er_cluster = calc_er_cluster(wt = wt, er = c(A, C, D, E, `F`, G, H, I, K, L, M, N, P, Q, R, S, `T`, V, W, Y)),
         pca_cluster = calc_pca_cluster(wt = wt, pca = c(PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14,
                                                         PC15, PC16, PC17, PC18, PC19, PC20))) %>%
  separate(er_cluster, c('er_dist', 'er_cluster'), sep = '_') %>%
  separate(pca_cluster, c('pca_dist', 'pca_cluster'), sep = '_') %>%
  mutate(er_perm_cluster = ifelse(is.na(perm_cluster), er_cluster, perm_cluster),
         pca_perm_cluster = ifelse(is.na(perm_cluster), pca_cluster, perm_cluster)) %>%
  select(study, gene, position, wt, cluster, er_cluster, er_perm_cluster, pca_cluster, pca_perm_cluster, er_dist, pca_dist)

# Analyse which method performs best, what distance suggests bad assignment, which clusters assign best/worst
