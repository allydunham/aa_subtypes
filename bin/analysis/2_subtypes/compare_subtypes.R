#!/usr/bin/env Rscript
# Compare properties between various different clustering schemes
source('src/config.R')
source('src/subtypes.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Compare AA subtypes', name = 'AA Subtype Comparison')
parser <- add_argument(parser, arg = '--subtypes', help = 'TSV file(s) assigning positions to clusters', nargs = Inf)
parser <- add_argument(parser, arg = '--dms', help = 'Path to DMS data', default = 'data/combined_mutational_scans.tsv')
parser <- add_argument(parser, arg = '--figures', help = 'Directory to save figures', default = '.')
args <- parse_args(parser)

subtypes <- sapply(args$subtypes, read_tsv, simplify = FALSE) %>%
  bind_rows(.id = 'method') %>%
  mutate(method = file_path_sans_ext(basename(method))) %>%
  tidyr::extract(method, into = c('algorithm', 'columns'), 
                 regex = "(dbscan|gmm|hclust|hdbscan|kmeans|pam)_(pca_no_sig|pca|profile)[[:alnum:]]*",
                 remove = FALSE)

dms <- read_tsv(args$dms) %>%
  left_join(subtypes, ., by = c('study', 'gene', 'position', 'wt'))

plots <- list()

### Silhouette Scores ###
silhouette_scores <- group_by(dms, method, wt) %>%
  group_modify(~mutate(., silhouette_score = cluster_silhouette(., A:Y)))

plots$method_silhouettes <- (ggplot(silhouette_scores, aes(x = method, y = silhouette_score, fill = algorithm)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = 'Silhouette Score (within AA)', x = 'Clustering Method') +
  guides(fill = guide_legend(title = 'Algorithm')) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dotted', colour = 'grey'))) %>%
  labeled_plot(units='cm', height=20, width=20)

### Cosine Similarity ###
get_cosine_sim <- function(tbl){
  cluster_profs <- cluster_mean_profiles(tbl, A:Y) %>%
    tibble_to_matrix(-cluster, row_names = 'cluster')
  
  # Cosine Similarity
  combs <- combn(nrow(cluster_profs), 2)
  tibble(cluster1 = rownames(cluster_profs)[combs[1,]],
         cluster2 = rownames(cluster_profs)[combs[2,]],
         cosine_sim = row_cosine_similarity(cluster_profs[combs[1,],], cluster_profs[combs[2,],])) %>%
    mutate(wt1 = str_sub(cluster1, end = 1),
           wt2 = str_sub(cluster2, end = 1)) %>%
    filter(wt1 == wt2)
}

cosine_sim <- filter(dms, !str_ends(cluster, '0')) %>%
  group_by(method) %>%
  group_modify(~get_cosine_sim(.)) %>%
  tidyr::extract(method, into = c('algorithm', 'columns'), 
                 regex = "(dbscan|gmm|hclust|hdbscan|kmeans)_(pca_no_sig|pca|profile)[[:alnum:]]*",
                 remove = FALSE)

plots$method_cosine_sim <- (ggplot(cosine_sim, aes(x = method, y = abs(cosine_sim), fill = algorithm)) +
  geom_boxplot() +
  geom_point() +
  coord_flip() +
  labs(y = 'Abs Cosine Similarity (within AA)', x = 'Clustering Method') +
  guides(fill = guide_legend(title = 'Algorithm')) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dotted', colour = 'grey'))) %>%
  labeled_plot(units='cm', height=20, width=20)

### Save plots ###
save_plotlist(plots, args$figures)
