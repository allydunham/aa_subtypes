#!/usr/bin/env Rscript
# Functions for AA subtypes clustering

#### k-means ####
make_kmeans_clusters <- function(tbl, cols, n=3, ...){
  cols <- enquo(cols)
  
  mat <- tibble_to_matrix(tbl, !!cols)
  
  km <- kmeans(mat, centers = n, ...)
  
  return(list(tbl=mutate(tbl, cluster = km$cluster) %>% select(cluster, everything()),
              kmeans=km))
}
########

#### Cluster Analysis ####
# All functions expect a cluster tibble
plot_ramachandran_angles <- function(tbl){
  tbl <- mutate(tbl, cluster_num = str_sub(cluster, start = -1))
  ggplot(tbl, aes(x=phi, y=psi, colour=cluster_num)) +
    geom_point() +
    facet_wrap(~wt) +
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(breaks = c(-180, -90, 0, 90, 180)) +
    scale_colour_viridis_d() +
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
    pivot_longer(-cluster, names_to = 'mut', values_to = 'score')
  
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
  cors <- cluster_profile_correlation(tbl, !!cols)
  
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

########