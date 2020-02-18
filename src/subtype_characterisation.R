#!/usr/bin/env Rscript
# Functions for AA subtypes characterisation
source('src/subtype_clustering.R')

#### Kmeans ####
evaluate_k_silhouette <- function(tbl, cols, k, dist_cols = NULL, min_size = 1){
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
  avg_sil <- filter(cluster_tbl, !str_detect(cluster, '^[A-Z]0$')) %>%
    group_by(wt) %>%
    group_modify(~mutate(., silhouette_score = cluster_silhouette(., !!dist_cols))) %>%
    summarise(silhouette_score = mean(silhouette_score, na.rm=TRUE))
  
  return(avg_sil)
}

evaluate_k_cosine <- function(tbl, cols, k, min_size = 1){
  cols <- enquo(cols)
  
  # Make Clusters
  clusters <- group_by(tbl, wt) %>%
    group_map(~make_kmeans_clusters(., !!cols, k = k, min_size = min_size), keep = TRUE) %>%
    set_names(sapply(., function(x){first(x$tbl$wt)}))
  
  cluster_profs <- map_dfr(clusters, .f = ~ .$tbl) %>%
    mutate(cluster = str_c(wt, cluster)) %>%
    arrange(study, position) %>%
    select(cluster, study, position, wt, A:Y) %>%
    cluster_mean_profiles(A:Y) %>%
    tibble_to_matrix(-cluster, row_names = 'cluster')
  
  # Cosine Similarity
  combs <- combn(nrow(cluster_profs), 2)
  cosine_sim <- tibble(cluster1 = rownames(cluster_profs)[combs[1,]],
                       cluster2 = rownames(cluster_profs)[combs[2,]],
                       cosine_sim = row_cosine_similarity(cluster_profs[combs[1,],], cluster_profs[combs[2,],])) %>%
    mutate(wt1 = str_sub(cluster1, end = 1),
           wt2 = str_sub(cluster2, end = 1)) %>%
    filter(wt1 == wt2)
  
  return(cosine_sim)
}

evaluate_k_sd <- function(tbl, cols, k, min_size = 1){
  cols <- enquo(cols)
  
  # Make Clusters
  clusters <- group_by(tbl, wt) %>%
    group_map(~make_kmeans_clusters(., !!cols, k = k, min_size = min_size), keep = TRUE) %>%
    set_names(sapply(., function(x){first(x$tbl$wt)}))
  
  cluster_profs <- map_dfr(clusters, .f = ~ .$tbl) %>%
    mutate(cluster = str_c(wt, cluster)) %>%
    arrange(study, position) %>%
    select(cluster, study, position, wt, A:Y) %>%
    cluster_mean_profiles(A:Y)
  
  # Cosine Similarity
  cluster_sd <- pivot_longer(cluster_profs, -cluster, names_to = 'aa', values_to = 'er') %>%
    group_by(cluster) %>%
    summarise(sd = sd(er)) %>%
    mutate(wt = str_sub(cluster, end = 1)) %>%
    group_by(wt) %>%
    summarise(sd = mean(sd))
  
  return(cluster_sd)
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

cluster_chem_env_profiles <- function(tbl, cols, normalise=NULL){
  cols <- enquo(cols)
  col_names <- colnames(select(tbl, !!cols))
  
  if (!is.null(normalise)){
    tbl <- mutate_at(tbl, vars(!!cols), as_mapper(normalise))
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
    mutate(rel_prob = log2(prob / mean(prob[!str_detect(cluster, '^[A-Z]0$')])))
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
  
  chem_env_profiles <- cluster_chem_env_profiles(tbl, within_10_0_A:within_10_0_Y, normalise = ~./max(., na.rm=TRUE)) %>%
    select(-n, -n_chem_env) %>%
    pivot_longer(within_10_0_A:within_10_0_Y, names_to = 'aa', names_prefix = 'within_10_0_', values_to = 'rel_count')
  
  aa_dist_normaliser <- function(x){
    x <- ifelse(is.finite(x), x, max(x[is.finite(x)]))
    x / max(x, na.rm = TRUE)
  }
  
  aa_distance_profiles <- cluster_chem_env_profiles(tbl, angstroms_to_A:angstroms_to_Y, normalise = aa_dist_normaliser) %>%
    select(-n, -n_chem_env) %>%
    pivot_longer(angstroms_to_A:angstroms_to_Y, names_to = 'aa', names_prefix = 'angstroms_to_', values_to = 'distance')
  
  sift_profiles <- cluster_sift_profiles(tbl) %>%
    pivot_longer(log10_sift_A:log10_sift_Y, names_to = 'aa', names_prefix = 'log10_sift_', values_to = 'log10_sift')
  
  ss_profiles <- cluster_ss_profile(tbl)
  
  sa_profiles <- cluster_sa_profile(tbl)
  
  l <- list(tbl=tbl, summary=cluster_summary, profiles=mean_profiles, sift=sift_profiles, foldx=foldx_profiles, 
            chem_env=chem_env_profiles, secondary_structure=ss_profiles, surface_accessibility=sa_profiles,
            aa_distance=aa_distance_profiles)
  class(l) <- c('cluster_characterisation', class(l))
  return(l)
}
########

#### Dendograms ####
# Add cluster labels to dendrogram d, based on position in c
add_clusters_to_dend <- function(d, c){
  dendrapply(d, function(x){attr(x, 'cluster') <- c[attr(x, 'label')]; x})
}

# Function to combine all nodes of the same cluster in a dendrogram
compress_dend <- function(d){
  # If we've recursed down to an individual leaf
  if (!is.null(attr(d, 'leaf'))){
    attr(d, 'label') <- str_c(attr(d, 'cluster'), ' (1)')
    return(d)
  }
  
  c <- squash(dendrapply(d, function(x){attr(x, 'cluster')}))
  
  # If all of the same cluster return a leaf
  if (all(map_lgl(c, ~. == c[[1]]))){
    leaf <- list()
    attributes(leaf) <- list(members=1, height=0, leaf=TRUE,
                             label=str_c(c[[1]], ' (', attr(d, 'members'), ')'),
                             cluster=c[[1]])
    return(leaf)
  } else {
    d[[1]] <- compress_dend(d[[1]])
    d[[2]] <- compress_dend(d[[2]])
    attr(d, 'members') <- attr(d[[1]], 'members') + attr(d[[2]], 'members')
    
    # Both leaves
    if (!is.null(attr(d[[1]], 'leaf')) & !is.null(attr(d[[2]], 'leaf'))){
      attr(d, 'midpoint') <- 0.5
      
      # Left leaf only
    } else if (!is.null(attr(d[[1]], 'leaf'))) {
      attr(d, 'midpoint') <- (1 + attr(d[[2]], 'midpoint'))/2
      
      # Right leaf only
    } else if (!is.null(attr(d[[2]], 'leaf'))) {
      attr(d, 'midpoint') <- attr(d[[1]], 'midpoint') + (attr(d[[1]], 'midpoint') + 1)/2
      
      # Both branches
    } else {
      attr(d, 'midpoint') <- attr(d[[1]], 'midpoint') + (attr(d[[1]], 'members') + attr(d[[2]], 'midpoint') + 1 - attr(d[[1]], 'midpoint'))/2
    }
    
    return(d)
  }
}

# Generic plot of a dendrogram from dendro_data type data (branches - tibble of branch data, leaves - tibble of leaf data with cluster assignment)
plot_dend <- function(branches, leaves){
  (ggplot() +
     geom_segment(data = branches, aes(x=x, y=y, xend=xend, yend=yend)) +
     geom_text(data = leaves, aes(x=x, y=y, label=label, colour=cluster), angle=90, hjust=1.2) +
     geom_point(data = leaves, aes(x=x, y=y, colour=cluster), shape=19) +
     scale_y_continuous(expand = expand_scale(mult = 0.15)) +
     theme(axis.line = element_blank(),
           axis.ticks = element_blank(),
           axis.text = element_blank(),
           axis.title = element_blank(),
           panel.grid.major.y = element_blank()) +
     guides(colour = guide_legend(title = 'Subtype', override.aes = list(label='', shape=15, size=3))) +
     scale_colour_brewer(type = 'qual', palette = 'Dark2', na.value = 'grey', direction = -1)) %>%
    labeled_plot(units='cm', height = 20, width = 20)
}

# Plot version of dendrograms with leaves combined together where a node is a single cluster
plot_compressed_dendrograms <- function(clusters, dms){
  dends <- map(names(clusters), ~add_clusters_to_dend(as.dendrogram(clusters[[.]]$hclust), filter(dms, wt == .) %>% pull(cluster))) %>%
    set_names(names(clusters))
  
  compressed_dends <- sapply(dends, compress_dend, simplify = FALSE)
  compressed_dend_data <- map(compressed_dends, dendro_data)
  
  branches <- map(compressed_dend_data,
                  ~as_tibble(.$segments) %>% 
                    mutate(yend = pmax(yend - 0.9 * min(y), 0),
                           y = pmax(y - 0.9 * min(y), 0)))
  
  leaves <- map(compressed_dend_data,
                ~as_tibble(.$labels) %>% 
                  mutate(label=as.character(label)) %>%
                  tidyr::extract(label, 'cluster', "([A-Z][0-9]*) \\([0-9]*\\)", remove = FALSE))
  
  map2(branches, leaves, plot_dend)
}

# Produce a dendogram from a given set of profiles, expected in wide format with cols cluster and A:Y
plot_profile_dendogram <- function(profiles, cols, distance_method='cosine'){
  cols <- enquo(cols)
  
  if (nrow(profiles) == 1){
    return(plot_dend(tibble(x=0, y=1, xend=0, yend=0),
                     tibble(x=0, y=0, label=profiles$cluster[1], cluster=profiles$cluster[1])))
  }
  
  if (distance_method == 'cosine'){
    distance <- tibble_to_matrix(profiles, !!cols, row_names = 'cluster') %>% cosine_distance_matrix() %>% as.dist()
  } else {
    distance <- tibble_to_matrix(profiles, !!cols, row_names = 'cluster') %>% dist(method = distance_method)
  }
  
  hc <- hclust(distance)
  dend_data <- dendro_data(hc)
  branches <- dend_data$segments
  leaves <- dend_data$labels
  
  multi_aa <- n_distinct(str_sub(leaves$label, end = 1)) > 1
  
  # Colour (cluster col in leaves) by AA if multiple AAs, otherwise by cluster
  if (multi_aa) {
    leaves <- mutate(leaves, cluster = str_sub(label, end = 1))
  } else {
    leaves <- mutate(leaves, cluster = label)
  }
  
  p <- plot_dend(branches, leaves)
  
  if (multi_aa){
    p <- p + scale_colour_manual(values = AA_COLOURS)
  }
  
  return(p)
}
########

#### Characterisation Plots ####
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
    x <- filter(x, !str_detect(cluster, '^[A-Z]0$'), n > filter_outliers)
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
    filter(!str_detect(cluster1, '^[A-Z]0$'),
           !str_detect(cluster2, '^[A-Z]0$')) %>%
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

plot_cluster_profile_distances <- function(x, method='manhattan'){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  distances <- cluster_mean_profiles(x, A:Y) %>%
    tibble_to_matrix(A:Y, row_names = 'cluster') %>%
    dist(method = method) %>%
    as.matrix() %>%
    as_tibble(rownames = 'cluster1') %>%
    pivot_longer(cols = -cluster1, names_to = 'cluster2', values_to = 'dist')
  
  ggplot(distances, aes(x = cluster1, y = cluster2, fill=dist)) +
    geom_tile() +
    scale_fill_distiller(type = ER_DIST_COLOURS$type, palette = ER_DIST_COLOURS$palette, direction = ER_DIST_COLOURS$direction) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = str_c(str_to_title(method), '\nDistance'))) +
    theme(axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(sort(unique(distances$cluster1)), end = 1)], angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(sort(unique(distances$cluster2)), end = 1)]))
}

plot_cluster_profile_cosine_sim <- function(x){
  if ('cluster_characterisation' %in% class(x)){
    x <- x$tbl
  }
  
  profs <- cluster_mean_profiles(x, A:Y) %>%
    tibble_to_matrix(A:Y, row_names = 'cluster')
  
  cosine_sim <- outer(1:nrow(profs), 1:nrow(profs), function(x, y){row_cosine_similarity(profs[x,], profs[y,])}) %>%
    set_colnames(rownames(profs)) %>%
    set_rownames(rownames(profs)) %>%
    as_tibble(rownames = 'cluster1') %>%
    pivot_longer(-cluster1, names_to = 'cluster2', values_to = 'cosine_sim') %>%
    filter(!str_detect(cluster1, '^[A-Z]0$'), !str_detect(cluster2, '^[A-Z]0$'))
  
  ggplot(cosine_sim, aes(x = cluster1, y = cluster2, fill=cosine_sim)) +
    geom_tile() +
    scale_fill_gradientn(colours = ER_COSINE_COLOURS$colours, limits = ER_COSINE_COLOURS$limits) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = 'Cosine\nSimilarity')) +
    theme(axis.ticks = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = AA_COLOURS[str_sub(sort(unique(cosine_sim$cluster1)), end = 1)], angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = AA_COLOURS[str_sub(sort(unique(cosine_sim$cluster2)), end = 1)]))
}

plot_cluster_foldx_profiles <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$foldx, select(x$summary, cluster, n, n_structure), by = 'cluster') %>%
      mutate(cluster = factor(cluster), prop = str_c(n_structure, '/', n))
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_detect(cluster, '^[A-Z]0$'), n_structure > filter_outliers)
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
    x <- filter(x, !str_detect(cluster, '^[A-Z]0$'), n_structure > filter_outliers)
  }
  
  cluster_labs <- levels(x$cluster)
  prop_labs <- structure(x$prop, names = as.character(x$cluster))[cluster_labs]
  
  breaks <- pretty_break(x$rel_count, rough_n = 5, sig_figs = 3)
  
  ggplot(x, aes(x=as.numeric(cluster), y=aa, fill=rel_count)) +
    geom_raster() +
    scale_fill_distiller(type = CHEM_ENV_COLOURS$type, palette = CHEM_ENV_COLOURS$palette, direction = CHEM_ENV_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    scale_x_continuous(breaks = 1:length(cluster_labs), labels = cluster_labs,
                       sec.axis = sec_axis(~., breaks = 1:length(prop_labs), labels = prop_labs)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression('log'[2]*frac('count', 'max(count)')))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(unique(x$aa))]),
          axis.text.x.bottom = element_text(colour = AA_COLOURS[str_sub(cluster_labs, end = 1)]),
          axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
          legend.title.align = 0.5)
}

plot_cluster_aa_distances <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$aa_distance, select(x$summary, cluster, n, n_structure), by = 'cluster') %>%
      mutate(cluster = factor(cluster), prop = str_c(n_structure, '/', n))
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_detect(cluster, '^[A-Z]0$'), n_structure > filter_outliers)
  }
  
  cluster_labs <- levels(x$cluster)
  prop_labs <- structure(x$prop, names = as.character(x$cluster))[cluster_labs]
  
  breaks <- pretty_break(x$distance, rough_n = 5, sig_figs = 3)
  
  ggplot(x, aes(x=as.numeric(cluster), y=aa, fill=distance)) +
    geom_tile() +
    scale_fill_distiller(type = AA_DISTANCE_COLOURS$type, palette = AA_DISTANCE_COLOURS$palette, direction = AA_DISTANCE_COLOURS$direction,
                         limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
    scale_x_continuous(breaks = 1:length(cluster_labs), labels = cluster_labs,
                       sec.axis = sec_axis(~., breaks = 1:length(prop_labs), labels = prop_labs)) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = expression(frac('||'~'Nearest X'~'||', '||'~'Nearest X'~'||'['max'])))) + 
    theme(plot.title = element_text(hjust = 0.5, size=8),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[sort(unique(x$aa))]),
          axis.text.x.bottom = element_text(colour = AA_COLOURS[str_sub(cluster_labs, end = 1)]),
          axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
          legend.title.align = 0.5)
}

plot_cluster_ss_profile <- function(x, filter_outliers=5){
  if ('cluster_characterisation' %in% class(x)){
    x <- left_join(x$secondary_structure, select(x$summary, cluster, n), by = 'cluster')
  }
  
  if (filter_outliers > 0){
    x <- filter(x, !str_detect(cluster, '^[A-Z]0$'), n > filter_outliers)
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
########

#### Full Characterisation Plot ####
# Generates a single combined plot for a set of clusters
plot_full_characterisation <- function(clusters, data, exclude_outliers=TRUE, global_scale=TRUE, outlier_size=10){
  # Clusters should all be sorted already in this workflow (1 largest etc.) but double check here as cheap to do,
  # plus have to add secondary tiebreaker anyway (so e.g. 9 and 10 are sorted as numbers not strings)
  cluster_order <- filter(data$summary, cluster %in% clusters) %>%
    mutate(cluster_num = as.integer(str_sub(cluster, start = 2))) %>%
    arrange(desc(n), cluster_num) %>%
    pull(cluster)
  
  # If only outliers exist just plot them, otherwise exclude clusters marked as outliers (X0) and under a given size
  if (exclude_outliers & !all(str_detect(cluster_order, '^[A-Z]0$'))){
    outliers <- filter(data$summary, cluster %in% clusters, str_detect(cluster, '^[A-Z]0$'), n < outlier_size) %>% pull(cluster)
    cluster_order <- cluster_order[!cluster_order %in% outliers]
  } else {
    outliers <- c()
  }
  
  # global outliers
  global_outliers <- filter(data$summary, str_detect(cluster, '^[A-Z]0$'), n < outlier_size) %>% pull(cluster)
  
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
  p_ss <- filter(data$secondary_structure, cluster %in% cluster_order, !ss=='i') %>%
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
  chem_env_lims <- filter(data$chem_env, cluster %in% cluster_order | global_scale, !cluster %in% global_outliers | !exclude_outliers) %>%  pull(rel_count) %>% pretty_break(rough_n = 4)
  p_chem_env <- filter(data$chem_env, cluster %in% cluster_order) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = cluster, y = aa, fill = rel_count)) +
    geom_raster() +
    coord_equal() +
    labs(y = '', x = '') +
    guides(fill = guide_colourbar(title = expression(frac('#AA', '#AA'[italic(max)])))) +
    scale_fill_distiller(type = CHEM_ENV_COLOURS$type, palette = CHEM_ENV_COLOURS$palette, direction = CHEM_ENV_COLOURS$direction,
                         limits = chem_env_lims$limits, breaks = chem_env_lims$breaks, labels = chem_env_lims$labels) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(colour = AA_COLOURS[unique(data$chem_env$aa)]),
          axis.text.x = element_text(colour = cluster_cols),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  text_theme <- theme(text = element_text(size = 9))
  legend_theme <- theme(legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'))
  
  if (length(cluster_order) < 7){
    size <- c(21, 21)
  } else {
    size <- c(29.7, 21)
  }
  
  p_overall <- multi_panel_figure(width = size[1], height = size[2], columns = 7, rows = 4, unit = 'cm', 
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
