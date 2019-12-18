#!/usr/bin/env Rscript
# Functions for AA subtypes clustering

#### Utility ####
# Lookup table for old style clustering scripts
CLUSTER_COLS <- list('profile'=quo(A:Y), 'pca'=quo(PC1:PC20), 'pca2'=quo(PC2:PC20))

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
  
  tbl <- mutate(tbl, cluster = compress_cluster_labels(cluster))
  
  return(list(tbl=tbl, kmeans=km))
}

# Expects named list of outputs from make_kmeans_clusters
plot_clustering_kmeans <- function(clusters){
  tbls <- sapply(clusters, extract2, 'tbl', simplify = FALSE) %>%
    bind_rows()
  
  ggplot(tbls, aes(x=umap1, y=umap2, colour=cluster)) +
    geom_point(shape = 20) +
    facet_wrap(~wt, nrow = 4) +
    labs(x='UMAP1', y='UMAP2') +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
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
  tbls <- sapply(clusters, extract2, 'tbl', simplify = FALSE) %>%
    bind_rows()
  
  ggplot(tbls, aes(x=umap1, y=umap2, colour=cluster)) +
    geom_point(shape = 20) +
    facet_wrap(~wt, nrow = 4) +
    labs(x='UMAP1', y='UMAP2') +
    guides(colour = guide_legend(title = 'Subtype')) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey')
}
########

#### Cluster Analysis ####
# Wrapper function
# Expects tbl to be in the format generated by make_dms_wide in subtypes_utils.R
# TODO - overall profile for cluster
make_cluster_plots <- function(tbl, cols, chem_env_cols, clusters){
  cols <- enquo(cols)
  chem_env_cols <- enquo(chem_env_cols)
  n_clusters <- n_distinct(tbl$cluster)
  plots <- list()
  
  plots$clustering <- labeled_plot(plot_clustering(clusters), units='cm', height = 40, width = 50)
  plots$tsne <- labeled_plot(plot_tsne_clusters(tbl), units='cm', height = 20, width = 20)
  plots$umap <- labeled_plot(plot_umap_clusters(tbl), units='cm', height = 20, width = 20)
  plots$ramachanran_angles <- labeled_plot(plot_ramachandran_angles(tbl), units='cm', height = 20, width = 20)
  plots$cluster_sizes <- labeled_plot(plot_cluster_sizes(tbl), units='cm', height = 10, width = n_clusters*0.5 + 2, limitsize=FALSE)
  plots$mean_profiles <- labeled_plot(plot_cluster_profiles(tbl, !!cols), units='cm', height = n_clusters*0.5 + 2, width = 15, limitsize=FALSE)
  plots$profile_correlation <- labeled_plot(plot_cluster_profile_correlation(tbl, !!cols), units='cm', height = n_clusters*0.5 + 2, width = n_clusters*0.5 + 4, limitsize=FALSE)
  plots$foldx_profiles <- labeled_plot(plot_cluster_foldx_profiles(tbl), units='cm', height = n_clusters*0.5 + 2, width = 25, limitsize=FALSE)
  plots$chem_env_profiles <- labeled_plot(plot_cluster_chem_env_profiles(tbl, !!chem_env_cols), units='cm', height = n_clusters*0.5 + 2, width = 25, limitsize=FALSE)
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
  } else {
    stop('Unrecognised clusters list')
  }
}

# All other functions expect a cluster tibble
plot_tsne_clusters <- function(tbl){
  mutate(dms_wide, cluster_sym = str_sub(cluster, start = -1)) %>%
    ggplot(aes(x=tSNE1, y=tSNE2, colour=cluster_sym)) +
    geom_point() +
    facet_wrap(~wt) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}

plot_umap_clusters <- function(tbl){
  mutate(dms_wide, cluster_sym = str_sub(cluster, start = -1)) %>%
    ggplot(aes(x=umap1, y=umap2, colour=cluster_sym)) +
    geom_point() +
    facet_wrap(~wt) +
    labs(x='UMAP1', y='UMAP2') +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    guides(colour = guide_legend(title = 'Subtype'))
}

plot_ramachandran_angles <- function(tbl){
  tbl <- mutate(tbl, cluster_num = str_sub(cluster, start = -1))
  ggplot(tbl, aes(x=phi, y=psi, colour=cluster_num)) +
    geom_point() +
    facet_wrap(~wt) +
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_colour_brewer(type = 'qual', palette = 'Set3', na.value = 'grey') +
    labs(x = expression(Phi), y = expression(Psi)) +
    guides(colour = guide_legend(title = 'Cluster')) +
    theme(panel.grid.major = element_line(linetype = 'dotted', colour='gray'))
}

plot_cluster_sizes <- function(tbl){
  size <- group_by(tbl, cluster) %>%
    tally() %>%
    mutate(aa = str_sub(cluster, end=1))
    
  ggplot(size, aes(x=cluster, y=n, fill=aa)) +
    geom_col() +
    scale_fill_manual(values = AA_COLOURS) +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(x='', y = 'Cluster Size')
}

cluster_mean_profiles <- function(tbl, cols){
  cols <- enquo(cols)
  group_by(tbl, cluster) %>%
    summarise_at(.vars = vars(A:Y), .funs = mean)
}

plot_cluster_profiles <- function(tbl, cols){
  cols <- enquo(cols)
  profiles <- cluster_mean_profiles(tbl, !!cols) %>%
    pivot_longer(-cluster, names_to = 'mut', values_to = 'score') %>%
    filter(str_ends(cluster, '0', negate = TRUE))
  
  ggplot(profiles, aes(x=mut, y=cluster, fill=score)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ggtitle('Cluster average profiles') +
    guides(fill=guide_colourbar(title = 'Normalised ER')) +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[unique(profiles$mut)]),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(unique(profiles$cluster), end = 1)]))
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

plot_cluster_profile_correlation <- function(tbl, cols){
  cols <- enquo(cols)
  cors <- cluster_profile_correlation(tbl, !!cols) %>%
    filter(str_ends(cluster1, '0', negate = TRUE),
           str_ends(cluster2, '0', negate = TRUE)) %>%
    mutate(cluster1 = droplevels(cluster1), cluster2 = droplevels(cluster2))
  
  ggplot(cors, aes(x=cluster1, y=cluster2, fill=cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle('Correlation between mean ER profiles') +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster1), end = 1)], angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(levels(cors$cluster2), end = 1)]))
}

cluster_foldx_profiles <- function(tbl){
    full_join(group_by(tbl, cluster) %>% summarise_at(vars(total_energy:energy_ionisation), mean, na.rm=TRUE),
              group_by(tbl, cluster) %>% summarise(n = n(), n_foldx = sum(!is.na(total_energy))),
              by = 'cluster')
}

plot_cluster_foldx_profiles <- function(tbl){
  profiles <- cluster_foldx_profiles(tbl) %>%
    filter(str_ends(cluster, '0', negate = TRUE)) %>%
    pivot_longer(cols = total_energy:energy_ionisation, names_to = 'term', values_to = 'ddg') %>%
    group_by(term) %>%
    mutate(rel_ddg = ddg / max(abs(ddg), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(term = FOLDX_TERMS[term])
  
  # Group terms based on complete profiles only
  term_order <- group_by(profiles, cluster) %>%
    filter(all(!is.nan(rel_ddg))) %>%
    ungroup() %>%
    add_factor_order(cluster, term, rel_ddg) %>%
    pull(term) %>%
    levels()
  
  profiles <- mutate(profiles,
                     term = factor(term, levels = term_order),
                     cluster = factor(cluster),
                     prop = str_c(n_foldx, '/', n))
  
  cluster_labs <- levels(profiles$cluster)
  prop_labs <- structure(profiles$prop, names = as.character(profiles$cluster))[cluster_labs]
  
  ggplot(profiles,
         aes(x=term, y=as.numeric(cluster), fill=rel_ddg)) +
    geom_raster() +
    scale_fill_gradient2() +
    scale_y_continuous(breaks = 1:length(cluster_labs),
                       labels = cluster_labs,
                       sec.axis = sec_axis(~.,
                                           breaks = 1:length(prop_labs),
                                           labels = prop_labs)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression(frac(Delta*Delta*'G', 'max'['term']*'(|'*Delta*Delta*'G|)')))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y.left = element_text(colour = AA_COLOURS[str_sub(unique(profiles$cluster), end = 1)]),
          legend.title.align = 0.5)
}

cluster_chem_env_profiles <- function(tbl, cols){
  cols <- enquo(cols)
  col_names <- colnames(select(tbl, !!cols))
  full_join(group_by(tbl, cluster) %>% summarise_at(vars(!!cols), mean, na.rm=TRUE),
            group_by(tbl, cluster) %>% summarise(n = n(), n_chem_env = sum(!is.na(.data[[col_names[1]]]))),
            by = 'cluster')
}

# Currently assumes profile columns correspond to AA counts, would need to be updated with a completely new chem env scheme
plot_cluster_chem_env_profiles <- function(tbl, cols){
  cols <- enquo(cols)
  profiles <- cluster_chem_env_profiles(tbl, !!cols) %>%
    pivot_longer(!!cols, names_to = 'aa', values_to = 'count') %>%
    group_by(aa) %>%
    mutate(rel_count = log2((count + min(count[count > 0], na.rm = TRUE)/10)/mean(count, na.rm=TRUE))) %>%
    ungroup() %>%
    mutate(prop = str_c(n_chem_env, '/', n),
           aa = str_sub(aa, start = -1)) %>%
    filter(str_ends(cluster, '0', negate = TRUE))
  
  # Group terms based on complete profiles only
  aa_order <- group_by(profiles, cluster) %>%
    filter(all(!is.nan(rel_count))) %>%
    ungroup() %>%
    add_factor_order(cluster, aa, rel_count) %>%
    pull(aa) %>%
    levels()
  
  profiles <- mutate(profiles, cluster = factor(cluster), aa = factor(aa, levels = aa_order))
  
  cluster_labs <- levels(profiles$cluster)
  prop_labs <- structure(profiles$prop, names = as.character(profiles$cluster))[cluster_labs]
  
  scale_cols <- c('red', 'orange', 'white', 'blue', 'purple')
  max_val <- max(abs(profiles$rel_count), na.rm = TRUE)
  scale_points <- rescale01(c(-max_val, -max_val/2, 0, max_val/2, max_val))
  
  ggplot(profiles,
         aes(x=aa, y=as.numeric(cluster), fill=rel_count)) +
    geom_raster() +
    scale_fill_gradientn(colours = scale_cols, values = scale_points, limits = c(-max_val, max_val)) +
    scale_y_continuous(breaks = 1:length(cluster_labs),
                       labels = cluster_labs,
                       sec.axis = sec_axis(~.,
                                           breaks = 1:length(prop_labs),
                                           labels = prop_labs)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression('log'[2]*frac('count', 'mean(count)')))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = AA_COLOURS[levels(profiles$aa)]),
          axis.text.y.left = element_text(colour = AA_COLOURS[str_sub(unique(profiles$cluster), end = 1)]),
          legend.title.align = 0.5)
}

### Full characterisation
full_cluster_characterisation <- function(tbl){
  cluster_summary <- group_by(tbl, cluster) %>%
    summarise(n = n(),
              n_structure = sum(!is.na(total_energy)),
              mean_sift = mean(mean_sift),
              mean_er = mean(mean_score),
              mean_foldx = mean(total_energy, na.rm=TRUE))
  
  mean_profiles <- cluster_mean_profiles(tbl, A:Y) %>%
    pivot_longer(A:Y, names_to = 'mut', values_to = 'er')
  
  foldx_profiles <- cluster_foldx_profiles(tbl) %>% 
    select(-n, -n_foldx) %>%
    filter(str_ends(cluster, '0', negate = TRUE)) %>%
    pivot_longer(cols = total_energy:energy_ionisation, names_to = 'term', values_to = 'ddg') %>%
    group_by(term) %>%
    mutate(rel_ddg = ddg / max(abs(ddg), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(term = FOLDX_TERMS[term]) %>%
    add_factor_order(cluster, term, rel_ddg) %>%
    mutate(cluster = as.character(cluster))
  
  chem_env_profiles <- cluster_chem_env_profiles(tbl, within_10_0_A:within_10_0_Y) %>%
    select(-n, -n_chem_env) %>%
    filter(str_ends(cluster, '0', negate = TRUE)) %>%
    pivot_longer(within_10_0_A:within_10_0_Y, names_to = 'aa', names_prefix = 'within_10_0_', values_to = 'count') %>%
    group_by(aa) %>%
    mutate(rel_count = log2((count + min(count[count > 0], na.rm = TRUE)/10)/mean(count, na.rm=TRUE))) %>%
    ungroup()
    
  sift_profiles <- group_by(tbl, cluster) %>%
    summarise_at(vars(log10_sift_A:log10_sift_Y), mean, na.rm=TRUE) %>%
    filter(str_ends(cluster, '0', negate = TRUE)) %>%
    pivot_longer(log10_sift_A:log10_sift_Y, names_to = 'aa', names_prefix = 'log10_sift_', values_to = 'log10_sift')
  
  ss_profiles <- group_by(tbl, cluster) %>%
    summarise_at(vars(starts_with('ss_')), mean, na.rm=TRUE) %>%
    pivot_longer(-cluster, names_to = 'ss', names_prefix = 'ss_', values_to = 'prob') %>%
    add_factor_order(cluster, ss, prob) %>%
    mutate(cluster = as.character(cluster)) %>%
    group_by(ss) %>%
    mutate(rel_prob = log2(prob / mean(prob)))
  
  sa_profiles <- group_by(tbl, cluster) %>%
    summarise_at(vars(all_atom_abs:polar_rel), mean, na.rm=TRUE)
  
  list(tbl=tbl, summary=cluster_summary, profiles=mean_profiles, sift=sift_profiles, foldx=foldx_profiles, 
       chem_env=chem_env_profiles, secondary_structure=ss_profiles, surface_accessibility=sa_profiles)
}

plot_full_characterisation <- function(clusters, data, exclude_outliers=TRUE, global_scale=TRUE){
  cluster_order <- filter(data$summary, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    arrange(desc(n)) %>%
    pull(cluster)
  
  cluster_cols <- AA_COLOURS[str_sub(cluster_order, end=1)]
  
  er_lims <- filter(data$profiles, !str_ends(cluster, '0') | !exclude_outliers, cluster %in% clusters | global_scale) %>%   pull(er) %>% pretty_break(step = 0.5, sym = 0)
  foldx_lims <- filter(data$foldx, !str_ends(cluster, '0') | !exclude_outliers, cluster %in% clusters | global_scale) %>%  pull(rel_ddg) %>%  pretty_break(step = 1, sym = 0)
  chem_env_lims <- filter(data$chem_env, !str_ends(cluster, '0') | !exclude_outliers, cluster %in% clusters | global_scale) %>%  pull(rel_count) %>% pretty_break(step = 1, sym = 0)
  ss_lims <- filter(data$secondary_structure, !str_ends(cluster, '0') | !exclude_outliers, cluster %in% clusters | global_scale) %>%  pull(rel_prob) %>% pretty_break(step = 1, sym = 0)
  sift_lims <- filter(data$sift, !str_ends(cluster, '0') | !exclude_outliers, cluster %in% clusters | global_scale) %>%  pull(log10_sift) %>% pretty_break(step = 1)
  
  # Summarise subtype sizes
  p_sizes <- filter(data$summary, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    mutate(cluster = factor(cluster, levels = rev(cluster_order)),
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
  p_sa <- filter(data$tbl, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
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
  p_ss <- filter(data$secondary_structure, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x=cluster, y=ss, fill=rel_prob)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    scale_fill_distiller(type = 'div', palette = 'BrBG', direction = 1, limits = ss_lims$limits, breaks = ss_lims$breaks, labels = ss_lims$labels) +
    scale_y_discrete(labels = DSSP_CLASSES) + 
    guides(fill = guide_colourbar(title = expression('log'[2]~frac(italic(P)(SS~"|"~X), italic(P)(SS))))) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean ER profiles
  p_profile <- filter(data$profiles, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = mut, fill = er)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = 'Norm. ER')) + 
    scale_fill_distiller(type = 'div', palette = 'RdBu', direction = 1, limits = er_lims$limits, breaks = er_lims$breaks, labels = er_lims$labels) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$profiles$mut)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean sift profiles
  p_sift <- filter(data$sift, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = aa, fill = log10_sift)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    scale_fill_distiller(type = 'seq', palette = 'PuRd', direction = -1, limits = sift_lims$limits, breaks = sift_lims$breaks, labels = sift_lims$labels) +
    guides(fill = guide_colourbar(title = expression(log[10]*'SIFT'))) + 
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$profiles$mut)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # Subtype mean foldx profiles
  p_foldx <- filter(data$foldx, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers, !term == 'Total Energy') %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = term, fill = rel_ddg)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = expression(frac(Delta*Delta*G[italic(term)], "|"*Delta*Delta*G[italic(term)]^{italic(max)}*"|")))) + 
    scale_fill_distiller(type = 'div', palette = 'PiYG', direction = 1, limits = foldx_lims$limits, breaks = foldx_lims$breaks, labels = foldx_lims$labels) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(colour = cluster_cols),
          axis.text.y = element_text(colour = 'black'), # TODO colour code foldx terms?
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  # subtype mean chemical environment profiles
  p_chem_env <- filter(data$chem_env, cluster %in% clusters, !str_ends(cluster, '0') | !exclude_outliers) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = aa, fill = rel_count)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = expression("log"[2]~frac('#AA', symbol("\341")~'#AA'~symbol("\361"))))) +
    scale_fill_distiller(type = 'div', palette = 'PRGn', direction = 1, limits = chem_env_lims$limits, breaks = chem_env_lims$breaks, labels = chem_env_lims$labels) +
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
    fill_panel(p_chem_env + text_theme + legend_theme, row = c(3, 4), column = c(6, 7), )
  
  return(list(overall=p_overall, sizes=p_sizes, secondary_structure=p_ss, profile=p_profile,
              sift=p_sift, foldx=p_foldx, chemial_environment=p_chem_env))
}
########