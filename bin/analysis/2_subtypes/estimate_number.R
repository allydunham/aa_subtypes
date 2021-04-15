#!/usr/bin/env Rscript
# Estimate number of subtypes by reclustering with varying numbers of positions
source('src/config.R')

dms <- read_tsv("data/combined_mutational_scans.tsv") %>%
  select(study, gene, position, wt, A:Y, PC1:PC20)

deep_split <- c("A" =  0, "C" =  0, "D" =  1, "E" =  0, "F" =  1, "G" =  1, "H" =  0, "I" =  0, "K" =  1,
                "L" =  0, "M" =  0, "N" =  0, "P" =  1, "Q" =  1, "R" =  1, "S" =  0, "T" =  1, "V" =  0,
                "W" =  0, "Y" =  1)

cosine_distance_matrix <- function(x) {
  return(acos(round(tcrossprod(x) / sqrt(tcrossprod(rowSums(x^2))), digits = 8)) / pi)
}

make_clusters <- function(tbl, deepSplit){
  mat <- tibble_to_matrix(tbl, PC2:PC20)
  d <- as.dist(cosine_distance_matrix(mat))
  hc <- hclust(d, method = "average")
  clus <- cutreeHybrid(dendro=hc, distM=as.matrix(d), deepSplit = deepSplit, verbose = 0)
  tbl <- mutate(tbl, cluster = as.character(clus$labels)) %>% select(cluster, everything())
  return(tbl)
}

subset_subtype_count <- function(tbl) {
  permissive_positions <- tibble_to_matrix(tbl, A:Y) %>%
    abs() %>%
    is_less_than(0.4) %>%
    apply(1, all)
  
  tbl_perm <- filter(tbl, permissive_positions) %>%
    mutate(cluster = CLUSTER_PERMISSIVE_CHAR)
  
  clusters <- filter(tbl, !permissive_positions) %>%
    group_by(wt) %>% 
    group_modify(~make_clusters(tbl = .x, deepSplit = deep_split[.y$wt])) %>%
    bind_rows(tbl_perm)
  
  counts <- filter(clusters, !cluster == "0") %>% 
    group_by(wt) %>%
    summarise(permissive = CLUSTER_PERMISSIVE_CHAR %in% cluster,
              nclusters = n_distinct(cluster[!cluster == CLUSTER_PERMISSIVE_CHAR]))
  return(counts)
}

increment_clustering <- function(start = 1000, inc = 200) {
  tot <- nrow(dms)
  shuf <- dms[sample(tot),]
  counts <- tibble(positions = unique(c(seq(start, tot, by = inc), tot))) %>%
    group_by(positions) %>%
    group_modify(~subset_subtype_count(shuf[seq_len(.y$positions),]))
  return(counts)
}

counts <- replicate(100, increment_clustering(), simplify = FALSE) %>%
  bind_rows(.id = "rep")
write_tsv(counts, "data/subtypes/rep_counts.tsv")

# Model number of subtypes as a logarithmic
models <- group_by(counts, wt) %>%
  {
    k <- group_keys(.)$wt
    group_map(., ~lm(nclusters ~ positions, data = .x)) %>%
      set_names(k)
  }
  

# Plot
plots <- list()
plots$subtype_count_boxes <- (ggplot(counts, aes(x = positions, group = positions, y = nclusters)) +
  facet_wrap(~wt) +
  geom_boxplot() +
  scale_y_continuous(breaks = 0:12) +
  labs(x = "Total number of positions (all amino acids)", y = "Functional Subtypes")) %>%
  labeled_plot(units = "cm", height = 25, width = 25)

plots$subtype_count_points <- (ggplot(counts, aes(x = positions, group = positions, y = nclusters)) +
  facet_wrap(~wt) +
  geom_jitter() +
  scale_y_continuous(breaks = 0:12) +
  labs(x = "Total number of positions (all amino acids)", y = "Functional Subtypes")) %>%
  labeled_plot(units = "cm", height = 25, width = 25)

save_plotlist(plots, "figures/2_subtypes/", overwrite = "all")
