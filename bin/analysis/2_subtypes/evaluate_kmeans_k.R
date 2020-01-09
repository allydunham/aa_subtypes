#!/usr/bin/env Rscript
# Determine optimal k for k-means clustering
source('src/config.R')
source('src/subtypes.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

k <- 2:7
avg_sil <- bind_rows(`Profile` = map(k, ~evaluate_k(dms, A:Y, ., min_size = 5)) %>% bind_rows(.id = 'k'),
                     `PCA` = map(k, ~evaluate_k(dms, PC1:PC20, ., min_size = 5)) %>% bind_rows(.id = 'k'),
                     `PCA (excluding PC1)` = map(k, ~evaluate_k(dms, PC2:PC20, ., min_size = 5)) %>% bind_rows(.id = 'k'),
                     .id = 'cols') %>%
  mutate(k = as.integer(k) + 1)

p_avg_sil <- ggplot(avg_sil, aes(x = k, y = silhouette_score, colour = cols)) +
  facet_wrap(~wt, scales = 'free', ncol = 5) +
  geom_point() +
  geom_line() +
  scale_color_brewer(type = 'qual', palette = 'Dark2', guide = guide_legend(title = '')) +
  labs(x = 'k', y = 'Silhouette Score (within AA)')
ggsave('figures/2_subtypes/kmean_k_silhouettes.pdf', p_avg_sil, units = 'cm', height = 20, width = 25)
