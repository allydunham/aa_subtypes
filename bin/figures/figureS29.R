#!/usr/bin/env Rscript
# Figure S29 (SIFT subtypes)
source('src/config.R')

dms <- left_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                  read_tsv('data/subtypes/sift_scores.tsv') %>% rename(sift_cluster = cluster),
                  by = c('study', 'gene', 'position', 'wt')) %>%
  left_join(read_tsv('data/combined_mutational_scans.tsv'), by = c('study', 'gene', 'position', 'wt'))

dms_dist <- filter(dms, !str_detect(cluster, CLUSTER_OUTLIER_RE)) %>%
  group_by(wt, cluster) %>%
  summarise_at(vars(A:Y), mean) %>%
  group_map(~tibble_to_matrix(., A:Y) %>% 
              cosine_distance_matrix() %>% 
              rowMeans() %>% 
              set_names(.x$cluster)) %>%
  unlist()

sift_dist <- filter(dms, !str_detect(sift_cluster, CLUSTER_OUTLIER_RE)) %>%
  group_by(wt, sift_cluster) %>%
  summarise_at(vars(starts_with('log10_sift')), mean) %>%
  group_map(~tibble_to_matrix(., starts_with('log10_sift')) %>% 
              cosine_distance_matrix() %>% 
              rowMeans() %>% 
              set_names(.x$sift_cluster)) %>%
  unlist()

cosine <- tibble(type=c(rep('dms', length(dms_dist)), rep('sift', length(sift_dist))),
                 cluster=c(names(dms_dist), names(sift_dist)),
                 cosine=c(dms_dist, sift_dist) %>% unname())

figure <- ggplot(cosine, aes(x = type, y = cosine, fill = type)) +
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(paired = FALSE, comparisons = list(c('dms', 'sift'))) +
  scale_x_discrete(name = 'Profiles used for clustering', labels = c(dms='ER', sift='log<sub>10</sub>SIFT4G')) +
  scale_y_continuous(name = 'Mean Cosine Distance', limits = c(0, 0.6)) +
  scale_fill_manual(values = c(dms='firebrick2', sift='cornflowerblue')) +
  theme(axis.text.x = element_markdown(),
        axis.ticks.x = element_blank())

ggsave('figures/4_figures/figureS29.pdf', figure, width = 120, height = 120, units = 'mm')
ggsave('figures/4_figures/figureS29.png', figure, width = 120, height = 120, units = 'mm')
ggsave('figures/4_figures/figureS29.tiff', figure, width = 120, height = 120, units = 'mm')
ggsave('figures/4_figures/figureS29.eps', figure, width = 120, height = 120, units = 'mm', device=cairo_ps, fallback_resolution = 600)

