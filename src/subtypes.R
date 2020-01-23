#!/usr/bin/env Rscript
# Functions for AA subtypes clustering
# TODO break down into smaller modules

#### Utility ####
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

# Calculate silhouette data for a clustering
# expects a tibble with columns cluster and cols, with cols giving the position
cluster_silhouette <- function(tbl, cols, distance_method = 'manhattan'){
  cols <- enquo(cols)
  
  if (n_distinct(tbl$cluster) == 1){
    warning('Silhouette undefined for a single cluster, returning NA')
    return(rep(NA, nrow(tbl)))
  }
  
  # Calculate distance between all points
  distance <- tibble_to_matrix(tbl, !!cols) %>%
    dist(method = distance_method) %>%
    as.matrix()
  diag(distance) <- NA
  
  # Calculate mean distance each point and all clusters
  get_silhouette_distance <- function(x){
    ind <- tbl$cluster == x
    if (sum(ind) == 1){
      return(distance[,ind])
    } else {
      return(rowMeans(distance[,tbl$cluster == x], na.rm = TRUE))
    }
  }
  mean_dists <- sapply(unique(tbl$cluster), get_silhouette_distance)
  
  # Calculate mean dist within cluster
  same_cluster <- cbind(1:nrow(mean_dists), match(tbl$cluster, colnames(mean_dists)))
  a <- mean_dists[same_cluster]
  
  # Calculate mean dist to other clusters
  mean_dists[same_cluster] <- NA
  b <- apply(mean_dists, 1, min, na.rm=TRUE)
  
  # Calculate silhouette score
  s <- (b - a)/(pmax(a, b))
  
  singlton_clusters <- table(tbl$cluster)
  singlton_clusters <- names(singlton_clusters)[singlton_clusters == 1]
  
  s[tbl$cluster %in% singlton_clusters] <- 0
  
  return(s)
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
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster))
  
  return(list(tbl=tbl, kmeans=km))
}

# Expects named list of outputs from make_kmeans_clusters
plot_clustering_kmeans <- function(clusters){
  tbls <- map_dfr(clusters, .f = ~ .$tbl) %>%
    bind_rows()
  
  ggplot(tbls, aes(x=umap1, y=umap2, colour=cluster)) +
    geom_point(shape = 20) +
    facet_wrap(~wt, nrow = 4) +
    labs(x='UMAP1', y='UMAP2') +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}

evaluate_k <- function(tbl, cols, k, dist_cols = NULL, min_size = 1){
  cols <- enquo(cols)
  dist_cols <- enquo(dist_cols)
  
  if (quo_is_null(dist_cols)){
    dist_cols <- cols 
  }
  
  # Make Clusters
  clusters <- group_by(tbl, wt) %>%
    group_map(~make_kmeans_clusters(., !!cols, k = k, min_size = min_size), keep = TRUE) %>%
    set_names(sapply(., function(x){first(x$tbl$wt)}))
  
  cluster_tbl <- map_dfr(clusters, .f = ~ .$tbl) %>%
    mutate(cluster = str_c(wt, cluster)) %>%
    arrange(study, position) %>%
    select(cluster, study, position, wt, !!dist_cols)
  
  # Average silhouette
  avg_sil <- filter(cluster_tbl, !str_ends(cluster, '0')) %>%
    group_by(wt) %>%
    group_modify(~mutate(., silhouette_score = cluster_silhouette(., !!dist_cols))) %>%
    summarise(silhouette_score = mean(silhouette_score, na.rm=TRUE))
    
  return(avg_sil)
}
########

#### hclust clustering ####
make_hclust_clusters <- function(tbl, cols, h = NULL, k = NULL, min_size = 1, dist_method = 'euclidean', method = 'average'){
  cols <- enquo(cols)

  mat <- tibble_to_matrix(tbl, !!cols)
  hc <- hclust(dist(mat, method = dist_method), method = method)
  clus <- cutree(hc, h = h, k = k)
  
  tbl <- mutate(tbl, cluster = clus) %>%
    select(cluster, everything())
  
  small_clusters <- count(tbl, cluster) %>%
    filter(n < min_size) %>%
    pull(cluster)
  
  tbl[tbl$cluster %in% small_clusters, 'cluster'] <- 0
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster))
  
  return(list(tbl = tbl, hclust = hc))
}

make_dynamic_hclust_clusters <- function(tbl, cols, dist_method = 'euclidean',
                                         hclust_args = list(method='average'),
                                         treecut_args = list()){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  d <- dist(mat, method = dist_method)
  hc <- do.call(hclust, c(list(d=d), hclust_args))
  clus <- do.call(cutreeHybrid, c(list(dendro=hc, distM=as.matrix(d)), treecut_args))
  
  tbl <- mutate(tbl, cluster = as.character(clus$labels)) %>%
    select(cluster, everything())
  
  return(list(tbl = tbl, hclust = hc))
}

# Expects named list of outputs from make_(dynamic_)hclust_clusters
plot_clustering_hclust <- function(clusters){
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
  
  tbl <- mutate(tbl, cluster = as.character(hdb$cluster)) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, hdbscan = hdb))
}

# Expects a named list of outputs from make_hdbscan_clusters
plot_clustering_hdbscan <- function(clusters){
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
  
  tbl <- mutate(tbl, cluster = as.character(db$cluster)) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, dbscan = db))
}
# Expects a named list of outputs from make_dbscan_clusters
plot_clustering_dbscan <- function(clusters){
  tbls <- map_dfr(clusters, .f = ~ .$tbl) %>%
    bind_rows()
  
  ggplot(tbls, aes(x=umap1, y=umap2, colour=cluster)) +
    geom_point(shape = 20) +
    facet_wrap(~wt, nrow = 4) +
    labs(x='UMAP1', y='UMAP2') +
    guides(colour = guide_legend(title = 'Subtype')) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey')
}
########

#### GMM clustering ####
make_gmm_clusters <- function(tbl, cols, G=1:5, modelNames = 'VVV', ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  gmm <- Mclust(mat, G, modelNames = modelNames)
  
  tbl <- mutate(tbl, cluster = as.character(gmm$classification)) %>% 
    select(cluster, everything())
  
  return(list(tbl = tbl, gmm = gmm))
}

# Expects a named list of outputs from make_gmm_clusters
get_mclustbic_matrix <- function(x){
  class(x) <- 'matrix'
  attributes(x) <- attributes(x)[c('dim', 'dimnames')]
  return(x)
}

plot_clustering_gmm <- function(clusters){
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

#### Characterise Clusters ####
# All characterisation functions expect a tibble in the format of make_dms_wide, with an additional cluster column
cluster_summary <- function(tbl){
  group_by(tbl, cluster) %>%
    summarise(n = n(),
              n_structure = sum(!is.na(total_energy)),
              mean_sift = mean(mean_sift),
              mean_er = mean(mean_score),
              mean_foldx = mean(total_energy, na.rm=TRUE),
              mean_sa = mean(all_atom_abs, na.rm=TRUE)) %>%
    mutate(aa = str_sub(cluster, end = 1))
}

cluster_mean_profiles <- function(tbl, cols){
  cols <- enquo(cols)
  group_by(tbl, cluster) %>%
    summarise_at(.vars = vars(A:Y), .funs = mean)
}

cluster_profile_correlation <- function(tbl, cols){
  cols <- enquo(cols)
  cluster_mean_profiles(tbl, !!cols) %>%
    transpose_tibble(cluster, id_col = 'aa') %>%
    tibble_correlation(x=-aa) %>%
    rename(cluster1 = cat1, cluster2 = cat2) %>%
    mutate(wt1 = str_sub(cluster1, end = 1),
           wt2 = str_sub(cluster2, end = 1))
}

cluster_sift_profiles <- function(tbl){
  group_by(tbl, cluster) %>%
    summarise_at(vars(log10_sift_A:log10_sift_Y), mean, na.rm=TRUE)
}

cluster_foldx_profiles <- function(tbl, normalise=FALSE){
  if (normalise){
    tbl <- mutate_at(tbl, vars(total_energy:energy_ionisation), ~./max(abs(.), na.rm = TRUE))
  }
  
  full_join(group_by(tbl, cluster) %>% summarise_at(vars(total_energy:energy_ionisation), mean, na.rm=TRUE),
            group_by(tbl, cluster) %>% summarise(n = n(), n_foldx = sum(!is.na(total_energy))),
            by = 'cluster')
}

cluster_chem_env_profiles <- function(tbl, cols, normalise=FALSE){
  cols <- enquo(cols)
  col_names <- colnames(select(tbl, !!cols))
  
  if (normalise){
    tbl <- mutate_at(tbl, vars(!!cols), ~log2((. + min(.[. > 0], na.rm = TRUE)/10)/mean(., na.rm=TRUE)), na.rm = TRUE)
  }
  
  full_join(group_by(tbl, cluster) %>% summarise_at(vars(!!cols), mean, na.rm=TRUE),
            group_by(tbl, cluster) %>% summarise(n = n(), n_chem_env = sum(!is.na(.data[[col_names[1]]]))),
            by = 'cluster')
}

# Gives average prob in a cluster, and avg. prob in a cluster relative to all non-noise clusters
cluster_ss_profile <- function(tbl){
  group_by(tbl, cluster) %>%
    summarise_at(vars(starts_with('ss_')), mean, na.rm=TRUE) %>%
    pivot_longer(-cluster, names_to = 'ss', names_prefix = 'ss_', values_to = 'prob') %>%
    add_factor_order(cluster, ss, prob) %>%
    mutate(cluster = as.character(cluster)) %>%
    group_by(ss) %>%
    mutate(rel_prob = log2(prob / mean(prob[!str_ends(cluster, '0')])))
}

cluster_sa_profile <- function(tbl){
  group_by(tbl, cluster) %>%
    summarise_at(vars(all_atom_abs:polar_rel), mean, na.rm=TRUE)
}

full_cluster_characterisation <- function(tbl){
  cluster_summary <- cluster_summary(tbl)
  
  mean_profiles <- cluster_mean_profiles(tbl, A:Y) %>%
    pivot_longer(A:Y, names_to = 'mut', values_to = 'er')
  
  foldx_profiles <- cluster_foldx_profiles(tbl, normalise = TRUE) %>% 
    pivot_longer(cols = total_energy:energy_ionisation, names_to = 'term', values_to = 'rel_ddg') %>%
    mutate(rel_ddg = ifelse(n_foldx == 0, 0, rel_ddg)) %>% # Have to set this to determine sensible term order, then unset
    add_factor_order(cluster, term, rel_ddg) %>%
    mutate(rel_ddg = ifelse(n_foldx == 0, NA, rel_ddg),
           cluster = as.character(cluster)) %>%
    select(-n, -n_foldx)
  
  chem_env_profiles <- cluster_chem_env_profiles(tbl, within_10_0_A:within_10_0_Y, normalise = TRUE) %>%
    select(-n, -n_chem_env) %>%
    pivot_longer(within_10_0_A:within_10_0_Y, names_to = 'aa', names_prefix = 'within_10_0_', values_to = 'rel_count')
  
  sift_profiles <- cluster_sift_profiles(tbl) %>%
    pivot_longer(log10_sift_A:log10_sift_Y, names_to = 'aa', names_prefix = 'log10_sift_', values_to = 'log10_sift')
  
  ss_profiles <- cluster_ss_profile(tbl)
  
  sa_profiles <- cluster_sa_profile(tbl)
  
  l <- list(tbl=tbl, summary=cluster_summary, profiles=mean_profiles, sift=sift_profiles, foldx=foldx_profiles, 
            chem_env=chem_env_profiles, secondary_structure=ss_profiles, surface_accessibility=sa_profiles)
  class(l) <- c('cluster_characterisation', class(l))
  return(l)
}
########

#### Plot/Analyse Cluster Characterisation ####
# Wrapper function
# Expects tbl to be in the format generated by make_dms_wide in subtypes_utils.R, with an additional clusters column
plot_cluster_diagnostics <- function(clusters, cols){
  cols <- enquo(cols)
  tbl <- map_dfr(clusters, .f = ~ .$tbl) %>%
    mutate(cluster = str_c(wt, cluster)) %>%
    arrange(study, position)
  
  n_clusters <- n_distinct(tbl$cluster)
  
  plots <- list()
  plots$clustering <- labeled_plot(plot_clustering(clusters), units='cm', height = 40, width = 50)
  plots$umap <- labeled_plot(plot_cluster_umap(tbl), units = 'cm', height = 20, width = 20)
  plots$tsne <- labeled_plot(plot_cluster_tsne(tbl), units = 'cm', height = 20, width = 20)
  plots$global_silhouette <- labeled_plot(plot_silhouette(tbl, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  plots$per_aa_silhouette <- labeled_plot(plot_per_aa_silhouette(tbl, A:Y), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  plots$clustering_silhouette <- labeled_plot(plot_per_aa_silhouette(tbl, !!cols), units='cm', height = n_clusters*0.33 + 2, width = 15, limitsize=FALSE)
  return(plots)
}

# Expects clusters as a list whose entries each have tbl and cluster items, as output by make_xx_clusters
plot_clustering <- function(clusters){
  # Dispatch to specific plot functions based on first entry of clusters list
  # (assume all will be the same as generated in this pipeline)
  if ('kmeans' %in% names(clusters[[1]])){
    p <- plot_clustering_kmeans(clusters)
  } else if ('hclust' %in% names(clusters[[1]])){
    p <- plot_clustering_hclust(clusters)
  } else if ('dbscan' %in% names(clusters[[1]])){
    p <- plot_clustering_dbscan(clusters)
  } else if ('hdbscan' %in% names(clusters[[1]])){
    p <- plot_clustering_hdbscan(clusters)
  } else if ('gmm' %in% names(clusters[[1]])){
    p <- plot_clustering_gmm(clusters)
  } else {
    stop('Unrecognised clusters list')
  }
}

# Silhouette plots work directly on a tibble, with a clusters ID column and cols defining spatial positions
# TODO refactor getting the silhouette score tbl out of these?
plot_silhouette <- function(tbl, cols){
  cols <- enquo(cols)
  
  tbl <- filter(tbl, !str_ends(cluster, '0')) %>%
    mutate(silhouette_score = cluster_silhouette(., !!cols)) %>%
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

plot_per_aa_silhouette <- function(tbl, cols){
  cols <- enquo(cols)
  
  tbl <- filter(tbl, !str_ends(cluster, '0')) %>%
    group_by(wt) %>%
    group_modify(~mutate(., silhouette_score = cluster_silhouette(., !!cols))) %>%
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

# Other functions utilise a part of a clustering full_characterisation, and can take either the subsection directly or full list
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

plot_cluster_ramachandran_angles <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  mutate(x, cluster_num = str_sub(cluster, start = -1)) %>%
  ggplot(aes(x=phi, y=psi, colour=cluster_num)) +
    geom_point() +
    facet_wrap(~wt) +
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_colour_brewer(type = CLUSTER_COLOURS$type, palette = CLUSTER_COLOURS$palette, direction = CLUSTER_COLOURS$direction, na.value = 'grey') +
    labs(x = expression(Phi), y = expression(Psi)) +
    guides(colour = guide_legend(title = 'Cluster')) +
    theme(panel.grid.major = element_line(linetype = 'dotted', colour='gray'))
}

plot_cluster_sizes <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$summary
  }
  
  ggplot(x, aes(x=cluster, y=n, fill=aa)) +
    geom_col() +
    geom_errorbar(aes(ymin=n_structure, ymax=n_structure), colour='white', width=0.4) +
    scale_fill_manual(values = AA_COLOURS) +
    guides(fill=FALSE) +
    labs(x='', y = 'Cluster Size')
}

plot_cluster_profiles <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$profiles, select(x$summary, cluster, n), by='cluster')
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_ends(cluster, '0'), n > filter_outliers)
  }
  
  breaks <- pretty_break(x$er, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(x, aes(x=cluster, y=mut, fill=er)) +
    geom_tile() +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    coord_fixed() +
    guides(fill=guide_colourbar(title = 'Norm. ER')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(unique(x$mut))]),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(sort(unique(x$cluster)), end = 1)]))
}

plot_cluster_profile_correlation <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  cors <- cluster_profile_correlation(x, A:Y) %>%
    filter(str_ends(cluster1, '0', negate = TRUE),
           str_ends(cluster2, '0', negate = TRUE)) %>%
    mutate(cluster1 = droplevels(cluster1), cluster2 = droplevels(cluster2))
  
  breaks <- pretty_break(cors$cor, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(cors, aes(x=cluster1, y=cluster2, fill=cor)) +
    geom_tile() +
    scale_fill_distiller(type = ER_COR_COLOURS$type, palette = ER_COR_COLOURS$palette, direction = ER_COR_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = 'Pearson\nCorrelation')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster1), end = 1)], angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster2), end = 1)]))
}

plot_cluster_foldx_profiles <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$foldx, select(x$summary, cluster, n, n_structure), by = 'cluster') %>%
      mutate(cluster = factor(cluster), prop = str_c(n_structure, '/', n))
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_ends(cluster, '0'), n_structure > filter_outliers)
  }
  
  cluster_labs <- levels(x$cluster)
  prop_labs <- structure(x$prop, names = as.character(x$cluster))[cluster_labs]
  
  breaks <- pretty_break(x$rel_ddg, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(x, aes(x=as.numeric(cluster), y=term, fill=rel_ddg)) +
    geom_raster() +
    scale_fill_distiller(type = FOLDX_COLOURS$type, palette = FOLDX_COLOURS$palette, direction = FOLDX_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    scale_x_continuous(breaks = 1:length(cluster_labs), labels = cluster_labs,
                       sec.axis = sec_axis(~., breaks = 1:length(prop_labs), labels = prop_labs)) +
    scale_y_discrete(labels = sapply(FOLDX_TERMS_PLOTMATH, function(x){parse(text = x)})) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression(frac(Delta*Delta*'G', 'max'['term']*'(|'*Delta*Delta*'G|)')))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x.bottom = element_text(colour = AA_COLOURS[str_sub(sort(unique(x$cluster)), end = 1)]),
          axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
          legend.title.align = 0.5)
}

# Currently assumes profile columns correspond to AA counts, would need to be updated with a completely new chem env scheme
plot_cluster_chem_env_profiles <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$chem_env, select(x$summary, cluster, n, n_structure), by = 'cluster') %>%
      mutate(cluster = factor(cluster), prop = str_c(n_structure, '/', n))
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_ends(cluster, '0'), n_structure > filter_outliers)
  }
  
  cluster_labs <- levels(x$cluster)
  prop_labs <- structure(x$prop, names = as.character(x$cluster))[cluster_labs]
  
  breaks <- pretty_break(x$rel_count, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(x, aes(x=as.numeric(cluster), y=aa, fill=rel_count)) +
    geom_raster() +
    scale_fill_distiller(type = CHEM_ENV_COLOURS$type, palette = CHEM_ENV_COLOURS$palette, direction = CHEM_ENV_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    scale_x_continuous(breaks = 1:length(cluster_labs), labels = cluster_labs,
                       sec.axis = sec_axis(~., breaks = 1:length(prop_labs), labels = prop_labs)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression('log'[2]*frac('count', 'mean(count)')))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(unique(x$aa))]),
          axis.text.x.bottom = element_text(colour = AA_COLOURS[str_sub(sort(unique(x$cluster)), end = 1)]),
          axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
          legend.title.align = 0.5)
}

plot_cluster_ss_profile <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$secondary_structure, select(x$summary, cluster, n), by = 'cluster')
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_ends(cluster, '0'), n > filter_outliers)
  }
  
  breaks <- pretty_break(x$rel_prob, rough_n = 5, sig_figs = 3, sym = 0)
  
  ggplot(x, aes(x=cluster, y=ss, fill=rel_prob)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    scale_fill_distiller(type = SS_COLOURS$type, palette = SS_COLOURS$palette, direction = SS_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    scale_y_discrete(labels = sapply(DSSP_CLASSES_PLOTMATH, function(x){parse(text = x)})) + 
    guides(fill = guide_colourbar(title = expression('log'[2]~frac(italic(P)(SS~"|"~X), italic(P)(SS))))) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(sort(unique(x$cluster)), end = 1)]),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
}

plot_cluster_ss_density <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  select(x, cluster, wt, starts_with('ss_')) %>%
    pivot_longer(starts_with('ss_'), names_to = 'ss', names_prefix = 'ss_', values_to = 'prob') %>%
    mutate(ss = DSSP_CLASSES_PLOTMATH[ss]) %>%
    ggplot(aes(x = prob, fill = wt)) +
    facet_grid(cols = vars(cluster), rows = vars(ss), scales = 'free_y', labeller = labeller(.rows = label_parsed)) +
    geom_density(alpha = 0.75, colour = NA) +
    scale_fill_manual(values = AA_COLOURS, guide = FALSE) +
    labs(x = expression(italic(P)(SS~"|"~"Cluster")), y = 'Density')
}

# Generates a single combined plot for a set of clusters
plot_full_characterisation <- function(clusters, data, exclude_outliers=TRUE, global_scale=TRUE, outlier_size=10){
  cluster_order <- filter(data$summary, cluster %in% clusters) %>%
    arrange(desc(n)) %>%
    pull(cluster)
  
  # If only outliers exist just plot them, otherwise exclude clusters marked as outliers (X0) and under a given size
  if (exclude_outliers & !all(str_ends(cluster_order, '0'))){
    outliers <- filter(data$summary, cluster %in% clusters, str_ends(cluster, '0'), n < outlier_size) %>% pull(cluster)
    cluster_order <- cluster_order[!cluster_order %in% outliers]
  } else {
    outliers <- c()
  }
  
  # global outliers
  global_outliers <- filter(data$summary, str_ends(cluster, '0'), n < outlier_size) %>% pull(cluster)
  
  cluster_cols <- AA_COLOURS[str_sub(cluster_order, end=1)]
  
  # Summarise subtype sizes (always include outliers in the cluster list here, for reference)
  p_sizes <- filter(data$summary, cluster %in% c(outliers, cluster_order)) %>%
    mutate(cluster = factor(cluster, levels = c(outliers, rev(cluster_order))),
           sum_str = str_c(n, ' (', n_structure, ')')) %>%
    ggplot(aes(x=cluster, fill=str_sub(cluster, end = 1))) +
    geom_col(aes(y=n), width = 0.5) +
    geom_errorbar(aes(ymin=n_structure, ymax=n_structure), colour='white', width=0.4) +
    geom_text(aes(y = max(n) + 10, label = sum_str), hjust=0, size=3) +
    scale_fill_manual(values = AA_COLOURS) +
    guides(fill=FALSE) +
    coord_flip(clip = 'off') +
    labs(x='', y='Count') +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
          axis.text.y = element_text(colour = cluster_cols),
          axis.ticks.y = element_blank())
  
  # Histograms of subtype surface accessibility
  p_sa <- filter(data$tbl, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = all_atom_abs, fill = wt)) +
    facet_wrap(~cluster, ncol = 1, strip.position = 'left') +
    geom_density(alpha = 0.75, colour = NA) +
    scale_fill_manual(values = AA_COLOURS) +
    labs(x = 'All Atom Abs.', y = '') +
    guides(fill = FALSE) +
    theme(strip.placement = 'outside',
          strip.text.y = element_text(angle = 180, colour = cluster_cols),
          panel.spacing = unit(0.05, 'npc'))
  
  # Subtype secondary structure probability increases vs background
  ss_lims <- filter(data$secondary_structure, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%  pull(rel_prob) %>% pretty_break(step = 1, sym = 0)
  p_ss <- filter(data$secondary_structure, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x=cluster, y=ss, fill=rel_prob)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    scale_fill_distiller(type = SS_COLOURS$type, palette = SS_COLOURS$palette, direction = SS_COLOURS$direction,
                         limits = ss_lims$limits, breaks = ss_lims$breaks, labels = ss_lims$labels) +
    scale_y_discrete(labels = sapply(DSSP_CLASSES_PLOTMATH, function(x){parse(text = x)})) + 
    guides(fill = guide_colourbar(title = expression('log'[2]~frac(italic(P)(SS~"|"~X), italic(P)(SS))))) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean ER profiles
  er_lims <- filter(data$profiles, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%   pull(er) %>% pretty_break(step = 0.5, sym = 0)
  p_profile <- filter(data$profiles, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = mut, fill = er)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = 'Norm. ER')) + 
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                         limits = er_lims$limits, breaks = er_lims$breaks, labels = er_lims$labels) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$profiles$mut)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean sift profiles
  sift_lims <- filter(data$sift, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%  pull(log10_sift) %>% pretty_break(step = 1)
  p_sift <- filter(data$sift, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = aa, fill = log10_sift)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    scale_fill_distiller(type = SIFT_COLOURS$type, palette = SIFT_COLOURS$palette, direction = SIFT_COLOURS$direction,
                         limits = sift_lims$limits, breaks = sift_lims$breaks, labels = sift_lims$labels) +
    guides(fill = guide_colourbar(title = expression(log[10]*'SIFT'))) + 
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$profiles$mut)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean foldx profiles
  foldx_lims <- filter(data$foldx, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%  pull(rel_ddg) %>%  pretty_break(step = 1, sym = 0)
  p_foldx <- filter(data$foldx, cluster %in% cluster_order, !term == 'Total Energy') %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = term, fill = rel_ddg)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = expression(frac(Delta*Delta*G[italic(term)], "|"*Delta*Delta*G[italic(term)]^{italic(max)}*"|")))) + 
    scale_fill_distiller(type = FOLDX_COLOURS$type, palette = FOLDX_COLOURS$palette, direction = FOLDX_COLOURS$direction,
                         limits = foldx_lims$limits, breaks = foldx_lims$breaks, labels = foldx_lims$labels) +
    scale_y_discrete(labels = sapply(FOLDX_TERMS_PLOTMATH, function(x){parse(text = x)})) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(colour = cluster_cols),
          axis.text.y = element_text(colour = 'black'), # TODO colour code foldx terms?
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # subtype mean chemical environment profiles
  chem_env_lims <- filter(data$chem_env, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%  pull(rel_count) %>% pretty_break(step = 1, sym = 0)
  p_chem_env <- filter(data$chem_env, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = aa, fill = rel_count)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = expression("log"[2]~frac('#AA', symbol("\341")~'#AA'~symbol("\361"))))) +
    scale_fill_distiller(type = CHEM_ENV_COLOURS$type, palette = CHEM_ENV_COLOURS$palette, direction = CHEM_ENV_COLOURS$direction,
                         limits = chem_env_lims$limits, breaks = chem_env_lims$breaks, labels = chem_env_lims$labels) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$chem_env$aa)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  text_theme <- theme(text = element_text(size = 9))
  legend_theme <- theme(legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'))
  
  p_overall <- multi_panel_figure(width = 21, height = 21, columns = 7, rows = 4, unit = 'cm', 
                                  row_spacing = 0.3, column_spacing = 0.3, panel_label_type = 'upper-alpha') %>%
    fill_panel(p_sizes + text_theme, row = 1, column = c(1, 2)) %>%
    fill_panel(p_ss + text_theme + legend_theme, row = 2, column = c(1, 2, 3)) %>%
    fill_panel(p_profile + text_theme + legend_theme, row = c(1, 2), column = c(4, 5)) %>%
    fill_panel(p_sift + text_theme + legend_theme, row = c(1, 2), column = c(6, 7)) %>%
    fill_panel(p_sa + text_theme, row = c(3, 4), column = c(1, 2)) %>%
    fill_panel(p_foldx + text_theme + legend_theme, row = c(3, 4), column = c(3, 4, 5)) %>%
    fill_panel(p_chem_env + text_theme + legend_theme, row = c(3, 4), column = c(6, 7))
  
  return(list(overall=p_overall, sizes=p_sizes, secondary_structure=p_ss, profile=p_profile,
              sift=p_sift, foldx=p_foldx, chemial_environment=p_chem_env))
}
########