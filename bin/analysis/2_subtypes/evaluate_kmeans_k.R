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
cosine_sim <- bind_rows(`Profile` = map(k, ~evaluate_k_cosine(dms, A:Y, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        `PCA` = map(k, ~evaluate_k_cosine(dms, PC1:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        `PCA (excluding PC1)` = map(k, ~evaluate_k_cosine(dms, PC2:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                        .id = 'cols') %>%
  mutate(k = as.integer(k))

# p_cosine_dist <- ggplot(cosine_sim, aes(x = k, y = cosine_sim, group = k)) + 
#   facet_grid(rows = vars(cols), cols = vars(wt1)) +
#   geom_boxplot() + 
#   geom_point()

avg_cosine <- group_by(cosine_sim, cols, k, wt=wt1) %>%
  summarise(abs_cosine_sim = mean(abs(cosine_sim)),
            cosine_sim = mean(cosine_sim))

p_avg_cosine <- ggplot(avg_cosine, aes(x = k, y = cosine_sim, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Set1', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Mean Cosine Similarity (within AA)')
ggsave('figures/2_subtypes/kmean_k_cosine.pdf', p_avg_cosine, units = 'cm', height = 20, width = 28)

p_avg_abs_cosine <- ggplot(avg_cosine, aes(x = k, y = abs_cosine_sim, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Set1', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Mean Absolute Cosine Similarity (within AA)')
ggsave('figures/2_subtypes/kmean_k_abs_cosine.pdf', p_avg_abs_cosine, units = 'cm', height = 20, width = 28)


### Cluster SD ###
avg_sd <- bind_rows(`Profile` = map(k, ~evaluate_k_sd(dms, A:Y, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                    `PCA` = map(k, ~evaluate_k_sd(dms, PC1:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                    `PCA (excluding PC1)` = map(k, ~evaluate_k_sd(dms, PC2:PC20, ., min_size = 5)) %>% set_names(k) %>% bind_rows(.id = 'k'),
                    .id = 'cols') %>%
  mutate(k = as.integer(k))

p_avg_sd <- ggplot(avg_sd, aes(x = k, y = sd, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Set1', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Profile Standard Deviation (within AA)')
ggsave('figures/2_subtypes/kmean_k_sd.pdf', p_avg_sd, units = 'cm', height = 20, width = 28)
