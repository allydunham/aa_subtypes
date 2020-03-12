#!/usr/bin/env Rscript
# Compare deepSplit values for dynamic tree cutting
source('src/config.R')

dms <- left_join(rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.tsv'), cluster_ds0 = cluster),
                 rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.tsv'), cluster_ds1 = cluster),
                 by = c("study", "gene", "position", "wt")) %>%
  left_join(read_tsv('data/combined_mutational_scans.tsv'), by = c("study", "gene", "position", "wt")) %>%
  select(cluster_ds0, cluster_ds1, everything())

# Get the mean correlation for profiles in a tibble
get_mean_cor <- function(dms, method='spearman'){
  f <- function(x, ...){
    tibble_to_matrix(x, A:Y) %>%
      t() %>%
      cor(method = method) %>%
      tril() %>%
      mean()
  }
  
  group_map(dms, f) %>%
    unlist() %>%
    tibble(cluster=group_keys(dms)[[1]], mean_cor=.)
}

mean_cors <- bind_rows(.id = 'deepSplit',
                       "0"=group_by(dms, cluster_ds0) %>% get_mean_cor(),
                       "1"=group_by(dms, cluster_ds1) %>% get_mean_cor()) %>%
  mutate(wt = str_sub(cluster, end = 1)) %>%
  filter(!str_detect(cluster, CLUSTER_OUTLIER_RE), !str_detect(cluster, CLUSTER_PERMISSIVE_RE))

p_mean_cors <- ggplot(mean_cors, aes(x = deepSplit, y = mean_cor, colour = wt)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~wt) +
  scale_color_manual(values = AA_COLOURS) +
  guides(colour=FALSE) +
  labs(y = expression('Cluster Mean Spearmans'~rho))
ggsave('figures/2_subtypes/hclust_dynamic_deep_split_mean_corelation.pdf', units = 'cm', width = 15, height = 15)