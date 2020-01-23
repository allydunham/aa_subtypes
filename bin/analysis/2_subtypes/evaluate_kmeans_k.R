#!/usr/bin/env Rscript
# Determine optimal k for k-means clustering
source('src/config.R')
source('src/subtypes.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

k <- 2:10

### Average Silhouette Score ###
avg_sil <- bind_rows(`Profile` = map(k, ~evaluate_k_silhouette(dms, A:Y, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                     `PCA` = map(k, ~evaluate_k_silhouette(dms, PC1:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                     `PCA (excluding PC1)` = map(k, ~evaluate_k_silhouette(dms, PC2:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                     `PCA (profile distance)` = map(k, ~evaluate_k_silhouette(dms, PC1:PC20, ., dist_cols = A:Y, min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                     `PCA (excluding PC1, profile distance)` = map(k, ~evaluate_k_silhouette(dms, PC2:PC20, ., dist_cols = A:Y, min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                     .id = 'cols') %>%
  mutate(k = as.integer(k))

p_avg_sil <- ggplot(avg_sil, aes(x = k, y = silhouette_score, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Set1', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Silhouette Score (within AA)')
ggsave('figures/2_subtypes/kmean_k_silhouettes.pdf', p_avg_sil, units = 'cm', height = 20, width = 28)

### Cluster Cosine ###
avg_cosine <- bind_rows(`Profile` = map(k, ~evaluate_k_cosine(dms, A:Y, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        `PCA` = map(k, ~evaluate_k_cosine(dms, PC1:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        `PCA (excluding PC1)` = map(k, ~evaluate_k_cosine(dms, PC2:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        .id = 'cols') %>%
  mutate(k = as.integer(k))

p_avg_cosine <- ggplot(avg_cosine, aes(x = k, y = cosine_sim, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Set1', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Cosine Similarity (within AA)')
ggsave('figures/2_subtypes/kmean_k_cosine.pdf', p_avg_cosine, units = 'cm', height = 20, width = 28)
