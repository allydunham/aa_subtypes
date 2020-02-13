#!/usr/bin/env Rscript
# Functions for AA subtypes clustering

#### Utility ####
# Reassign cluster numbers by size
order_cluster_labels <- function(x){
  counts <- table(x)
  counts <- counts[names(counts) != '0']
  hash <- structure(as.character(1:length(counts)), names=names(counts)[order(counts, decreasing = TRUE)])
  hash['0'] <- '0'
  
  unname(hash[as.character(x)])
}

# Drop unused integer cluster labels
compress_cluster_labels <- function(x){
  unq <- unique(x)
  if (0 %in% unq){
    hash <- structure(as.character(0:(length(unq)-1)), names=as.character(sort(unq)))
  } else {
    hash <- structure(as.character(1:length(unq)), names=as.character(sort(unq)))
  }
  
  unname(hash[as.character(x)])
}
########

#### k-means ####
make_kmeans_clusters <- function(tbl, cols, k=3, min_size=1, ...){
  cols <- enquo(cols)

  mat <- tibble_to_matrix(tbl, !!cols)
  
  km <- kmeans(mat, centers = k, ...)
  
  tbl <- mutate(tbl, cluster = km$cluster) %>%
    select(cluster, everything())
  
  small_clusters <- count(tbl, cluster) %>%
    filter(n < min_size) %>%
    pull(cluster)
  
  tbl[tbl$cluster %in% small_clusters, 'cluster'] <- 0
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster) %>% order_cluster_labels())
  
  return(list(tbl=tbl, kmeans=km))
}

make_multi_kmeans_clusters  <- function(tbl, cols, k=3, min_size=1, min_k=2, reps=100){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  km_reps <- replicate(reps, get_kmeans_rep(mat, k, min_size, min_k), simplify = FALSE)
  mean_abs_cosim <- map_dbl(km_reps, ~mean(.$cosim))
  
  km <- km_reps[[which.min(mean_abs_cosim)]]$km
  
  tbl <- mutate(tbl, cluster = km$cluster) %>%
    select(cluster, everything())
  
  small_clusters <- count(tbl, cluster) %>%
    filter(n < min_size) %>%
    pull(cluster)
  
  tbl[tbl$cluster %in% small_clusters, 'cluster'] <- 0
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster) %>% order_cluster_labels())
  
  return(list(tbl=tbl, kmeans=km))
}

get_kmeans_rep <- function(mat, k, min_size=1, min_k=2){
  km <- kmeans(mat, centers = k)
  
  sizes <- table(km$cluster)
  large_clusters <- as.integer(names(sizes[sizes > min_size]))
  
  if (length(large_clusters) < min_k){
    return(list(km=km, cosim=Inf))
  }
  
  com <- combn(large_clusters, 2)
  
  if (length(large_clusters) == 2){
    x <- km$centers[com[1,],]
    y <- km$centers[com[2,],]
    cosim <- abs(sum(x*y) / (sqrt(sum(x^2) * sum(y^2))))
  } else {
    cosim <- abs(row_cosine_similarity(km$centers[com[1,],], km$centers[com[2,],]))
  }
  
  return(list(km=km, cosim=cosim))
}
########

#### PAM ####
make_pam_clusters <- function(tbl, cols, k=3, min_size=1, distance_method='cosine', ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  if (distance_method == 'cosine'){
    combs <- combn(1:nrow(mat), 2)
    d <- acos(row_cosine_similarity(mat[combs[1,],], mat[combs[2,],])) / pi
  } else {
    d <- dist(mat, method = distance_method)
  }
  
  pm <- pam(d, k = k, diss=TRUE, ...)
  
  tbl <- mutate(tbl, cluster = pm$clustering) %>%
    select(cluster, everything())
  
  small_clusters <- count(tbl, cluster) %>%
    filter(n < min_size) %>%
    pull(cluster)
  
  tbl[tbl$cluster %in% small_clusters, 'cluster'] <- 0
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster) %>% order_cluster_labels())
  
  return(list(tbl=tbl, pam=pm))
}

plot_pam_medoids <- function(clusters){
  tbls <- map_dfr(clusters, .f = ~ .$tbl) %>%
    bind_rows()
  
  medoids <- map(clusters, ~select(.$tbl, cluster, wt, umap1, umap2)[.$pam$medoids,]) %>%
    bind_rows()
  
  ggplot(tbls, aes(x=umap1, y=umap2, colour=cluster)) +
    geom_point(shape = 20) +
    geom_point(data = medoids, aes(x=umap1, y=umap2, fill=cluster), colour='black', shape=23, size=3) +
    facet_wrap(~wt, nrow = 4) +
    labs(x='UMAP1', y='UMAP2') +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    scale_fill_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}
########

#### hclust clustering ####
make_hclust_clusters <- function(tbl, cols, h = NULL, k = NULL, min_size = 1, distance_method = 'euclidean', method = 'average'){
  cols <- enquo(cols)

  mat <- tibble_to_matrix(tbl, !!cols)
  
  if (distance_method == 'cosine'){
    d <- as.dist(cosine_distance_matrix(mat))
  } else {
    d <- dist(mat, method = distance_method)
  }
  
  hc <- hclust(d, method = method)
  clus <- cutree(hc, h = h, k = k)
  
  tbl <- mutate(tbl, cluster = clus) %>%
    select(cluster, everything())
  
  small_clusters <- count(tbl, cluster) %>%
    filter(n < min_size) %>%
    pull(cluster)
  
  tbl[tbl$cluster %in% small_clusters, 'cluster'] <- 0
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster) %>% order_cluster_labels())
  
  return(list(tbl = tbl, hclust = hc))
}

make_dynamic_hclust_clusters <- function(tbl, cols, distance_method = 'euclidean',
                                         hclust_args = list(method='average'),
                                         treecut_args = list()){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  if (distance_method == 'cosine'){
    d <- as.dist(cosine_distance_matrix(mat))
  } else {
    d <- dist(mat, method = distance_method)
  }
  
  hc <- do.call(hclust, c(list(d=d), hclust_args))
  clus <- do.call(cutreeHybrid, c(list(dendro=hc, distM=as.matrix(d)), treecut_args))
  
  tbl <- mutate(tbl, cluster = as.character(clus$labels)) %>% # Already has cluster labels ordered by size
    select(cluster, everything())
  
  return(list(tbl = tbl, hclust = hc))
}

# Expects named list of outputs from make_(dynamic_)hclust_clusters
plot_hclust_dend <- function(clusters){
  dend_data <- sapply(clusters, function(x){dendro_data(x$hclust)}, simplify = FALSE)
  
  branches <- sapply(dend_data, extract2, 'segments', simplify = FALSE) %>%
    bind_rows(.id = 'wt') %>%
    as_tibble()
  
  leaves <- sapply(dend_data, function(x){mutate(as_tibble(x$labels), label = as.character(label))}, simplify = FALSE) %>%
    map2(clusters, function(x, y){bind_cols(x, y$tbl[as.integer(x$label),'cluster'])}) %>%
    bind_rows(.id = 'wt')
  
  ggplot() +
    geom_segment(data = branches, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_point(data = leaves, aes(x=x, y=y, colour=cluster), shape = 20) +
    facet_wrap(~wt, nrow = 4, scales = 'free_x') +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank()) +
    guides(colour = guide_legend(title = 'Subtype')) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey')
}
########

#### hdbscan clustering ####
make_hdbscan_clusters <- function(tbl, cols, dist_method = 'euclidean', minPts=10, ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  dis <- dist(mat, method = dist_method)
  hdb <- hdbscan(mat, minPts = minPts, xdist = dis, ...)
  
  tbl <- mutate(tbl, cluster = as.character(hdb$cluster) %>% order_cluster_labels()) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, hdbscan = hdb))
}

# Expects a named list of outputs from make_hdbscan_clusters
plot_hdbscan_dend <- function(clusters){
  dend_data <- sapply(clusters, function(x){dendro_data(x$hdbscan$hc)}, simplify = FALSE)
  
  branches <- sapply(dend_data, extract2, 'segments', simplify = FALSE) %>%
    bind_rows(.id = 'wt') %>%
    as_tibble()
  
  leaves <- sapply(dend_data, function(x){mutate(as_tibble(x$labels), label = as.character(label))}, simplify = FALSE) %>%
    map2(clusters, function(x, y){bind_cols(x, y$tbl[as.integer(x$label),'cluster'])}) %>%
    bind_rows(.id = 'wt')
  
  ggplot() +
    geom_segment(data = branches, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_point(data = leaves, aes(x=x, y=y, colour=cluster), shape = 20) +
    facet_wrap(~wt, nrow = 4, scales = 'free_x') +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank()) +
    guides(colour = guide_legend(title = 'Subtype')) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey')
}
########

#### dbscan clustering ####
make_dbscan_clusters <- function(tbl, cols, eps, dist_method = 'euclidean', minPts=5, ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  dis <- dist(mat, method = dist_method)
  db <- dbscan(dis, eps=eps, minPts = minPts, ...)
  
  tbl <- mutate(tbl, cluster = as.character(db$cluster) %>% order_cluster_labels()) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, dbscan = db))
}
########

#### GMM clustering ####
make_gmm_clusters <- function(tbl, cols, G=1:5, modelNames = 'VVV', ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  gmm <- Mclust(mat, G, modelNames = modelNames)
  
  tbl <- mutate(tbl, cluster = as.character(gmm$classification) %>% order_cluster_labels()) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, gmm = gmm))
}

# Expects a named list of outputs from make_gmm_clusters
get_mclustbic_matrix <- function(x){
  class(x) <- 'matrix'
  attributes(x) <- attributes(x)[c('dim', 'dimnames')]
  return(x)
}

plot_gmm_bic <- function(clusters){
  bic <- map(clusters, .f = ~ as_tibble(get_mclustbic_matrix(.$gmm$BIC), rownames = 'n')) %>%
    bind_rows(.id = 'aa') %>%
    pivot_longer(c(-aa, -n), names_to = 'model', values_to = 'BIC')
  
  ggplot(bic, aes(x=n, y=BIC, fill=model)) +
    geom_col(position = 'dodge') +
    facet_wrap(~aa, nrow = 4) +
    labs(x='Number of Subtypes', y='BIC') +
    guides(fill = guide_legend(title = 'Model')) +
    scale_fill_brewer(type = 'qual', palette = 'Set1', na.value = 'grey')
}
########

#### Plot Cluster Diagnostics ####
# Wrapper function
# Expects tbl to be in the format generated by make_dms_wide in subtypes_utils.R, with an additional clusters column
plot_cluster_diagnostics <- function(clusters, cols){
  cols <- enquo(cols)
  tbl <- map_dfr(clusters, .f = ~ .$tbl) %>%
    mutate(cluster = str_c(wt, cluster)) %>%
    arrange(study, position)
  
  n_clusters <- n_distinct(tbl$cluster)
  
  plots <- list()
  plots$umap <- labeled_plot(plot_cluster_umap(tbl), units = 'cm', height = 20, width = 20)
  plots$tsne <- labeled_plot(plot_cluster_tsne(tbl), units = 'cm', height = 20, width = 20)
  plots$silhouette_global <- labeled_plot(plot_silhouette(tbl, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  plots$silhouette_per_aa <- labeled_plot(plot_per_aa_silhouette(tbl, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  plots$silhouette_per_aa_cosine <- labeled_plot(plot_per_aa_silhouette(tbl, A:Y, 'cosine'), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  plots$silhouette_clustering_vars <- labeled_plot(plot_per_aa_silhouette(tbl, !!cols), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  
  if ('hclust' %in% names(clusters[[1]])){
    plots$hclust_dend <- plot_hclust_dend(clusters)
  } else if ('hdbscan' %in% names(clusters[[1]])){
    plots$hdbscan_dend <- plot_hdbscan_dend(clusters)
  } else if ('gmm' %in% names(clusters[[1]])){
    plots$gmm_bic <- plot_gmm_bic(clusters)
  } else if ('pam' %in% names(clusters[[1]])){
    plots$medoids <- plot_pam_medoids(clusters)
  }
  
  return(plots)
}

# Silhouette plots work directly on a tibble, with a clusters ID column and cols defining spatial positions
plot_silhouette <- function(tbl, cols, distance_method='manhattan'){
  cols <- enquo(cols)
  
  tbl <- filter(tbl, !str_ends(cluster, '0')) %>%
    mutate(silhouette_score = cluster_silhouette(., !!cols, distance_method = distance_method)) %>%
    select(cluster, study, position, wt, silhouette_score)
  
  ggplot(tbl, aes(x = cluster, y = silhouette_score, fill=wt)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = AA_COLOURS, guide = FALSE) +
    labs(x = '', y = 'Silhouette Score (global)') +
    coord_flip() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(str_sub(unique(tbl$cluster), end = 1))]))
}

plot_per_aa_silhouette <- function(tbl, cols, distance_method='manhattan'){
  cols <- enquo(cols)
  
  tbl <- filter(tbl, !str_ends(cluster, '0')) %>%
    group_by(wt) %>%
    group_modify(~mutate(., silhouette_score = cluster_silhouette(., !!cols, distance_method = distance_method))) %>%
    select(cluster, study, position, wt, silhouette_score) %>%
    ungroup()
  
  ggplot(tbl, aes(x = cluster, y = silhouette_score, fill=wt)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = AA_COLOURS, guide = FALSE) +
    labs(x = '', y = 'Silhouette Score (within AA)') +
    coord_flip() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(str_sub(unique(tbl$cluster), end = 1))]))
}

plot_cluster_umap <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  mutate(x, cluster_sym = str_sub(cluster, start = -1)) %>%
    ggplot(aes(x=umap1, y=umap2, colour=cluster_sym)) +
    geom_point() +
    facet_wrap(~wt) +
    labs(x='UMAP1', y='UMAP2') +
    scale_colour_brewer(type = CLUSTER_COLOURS$type, palette = CLUSTER_COLOURS$palette, direction = CLUSTER_COLOURS$direction, na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}

plot_cluster_tsne <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  mutate(x, cluster_sym = str_sub(cluster, start = -1)) %>%
    ggplot(aes(x=tSNE1, y=tSNE2, colour=cluster_sym)) +
    geom_point() +
    facet_wrap(~wt) +
    scale_colour_brewer(type = CLUSTER_COLOURS$type, palette = CLUSTER_COLOURS$palette, direction = CLUSTER_COLOURS$direction, na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}
########
