#!/ussr/bin/env Rscript
# Perform a K-means clustering of AA subtypes

library(argparser)
source('src/config.R')
source('src/clustering.R')

### Parse Args and perform setup ###
parser <- arg_parser(description = 'Make and analyse AA subtypes using kmeans clustering', name = 'Kmeans AA Clustering')
parser <- add_argument(parser, arg = '--ncluster', help = 'Number of clusters to split each AA into', default = 3)
args <- parse_args(parser)

print(args$ncluster)

root_name <- str_c('kmeans_', args$ncluster)
root_fig_dir <- str_c('figures/2_clustering/', root_name)
dir.create(root_fig_dir, recursive = TRUE)
dir.create('data/clusterings')
plots <- list()

### Import data ###
dms <- read_tsv('data/combined_mutational_scans.tsv')

foldx_averages <- select(dms, study, position, wt, total_energy:entropy_complex) %>%
  select(-sloop_entropy, -mloop_entropy, -entropy_complex, -water_bridge) %>% # Drop terms that are unused in our structures
  drop_na(total_energy) %>%
  group_by(study, position, wt) %>%
  summarise_all(mean, na.rm=TRUE)

dms_wide <- filter(dms, mut %in% Biostrings::AA_STANDARD) %>%
  select(study, gene, position, wt, mut, imputed_score, log10_sift, psi, phi) %>%
  pivot_wider(names_from = mut, values_from = c(imputed_score, log10_sift)) %>%
  rename_at(vars(starts_with('imputed_score_')), ~str_sub(., start=-1)) %>%
  mutate(mean_score = rowMeans(select(., A:Y)),
         mean_sift = rowMeans(select(., log10_sift_A:log10_sift_Y))) %>%
  left_join(foldx_averages, by = c('study', 'position', 'wt'))

### Create Clusters ###
kmeans_cluster <- group_by(dms_wide, wt) %>%
  group_map(~make_kmeans_clusters(., A:Y, n = args$ncluster), keep = TRUE)

dms_wide <- map_dfr(kmeans_cluster, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, cluster)) %>%
  arrange(study, position)

n_clusters <- n_distinct(dms_wide$cluster)

### Save cluster tsv ###
write_tsv(select(dms_wide, cluster, study, gene, position, wt), str_c('data/clusterings/', root_name, '.tsv'))

### Analyse clusters ###
plots$ramachanran_angles <- labeled_plot(plot_ramachandran_angles(dms_wide), units='cm', height = 20, width = 20)
plots$cluster_sizes <- labeled_plot(plot_cluster_sizes(dms_wide), units='cm', height = 10, width = 14)
plots$mean_profiles <- labeled_plot(plot_cluster_profiles(dms_wide, A:Y), units='cm', height = n_clusters*0.5 + 2, width = 15)
plots$profile_correlation <- labeled_plot(plot_cluster_profile_correlation(dms_wide, A:Y), units='cm', height = n_clusters*0.5 + 2, width = n_clusters*0.5 + 4)

save_plotlist(plots, root = root_fig_dir)
