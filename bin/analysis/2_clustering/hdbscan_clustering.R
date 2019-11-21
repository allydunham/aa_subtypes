#!/ussr/bin/env Rscript
# Perform an HDBSCAN clustering of AA subtypes

### Parse Args and perform setup ###
library(argparser)
parser <- arg_parser(description = 'Make and analyse AA subtypes using kmeans clustering', name = 'Kmeans AA Clustering')
parser <- add_argument(parser, arg = '--min_size', help = 'Minimum allowed cluster size', default = 10)
parser <- add_argument(parser, arg = '--mode', help = 'Cluster using "profile" or "pca"', default = 'profile')
args <- parse_args(parser)

if (!args$mode %in% c('profile', 'pca')){
  stop('--mode must be one of "profile" or "pca"')
}

library(dbscan)
source('src/config.R')
source('src/clustering.R')
root_name <- str_c('hdbscan', args$mode, args$min_size, sep = '_')
root_fig_dir <- str_c('figures/2_clustering/', root_name)
dir.create(root_fig_dir, recursive = TRUE)
dir.create('data/clusterings')

### Import data ###
dms <- read_tsv('data/combined_mutational_scans.tsv')
dms_wide <- make_dms_wide(dms)
pca <- tibble_pca(dms_wide, A:Y)
dms_wide <- bind_cols(dms_wide, as_tibble(pca$x))

### Create Clusters ###
if (args$mode == 'profile'){
  cols <- quo(A:Y)
} else if (args$mode == 'pca'){
  cols <- quo(PC2:PC20)
}

hdbscan_cluster <- group_by(dms_wide, wt) %>%
  group_map(~make_hdbscan_clusters(., !!cols, minPts = args$min_size), keep = TRUE)

dms_wide <- map_dfr(hdbscan_cluster, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, ifelse(cluster == 0, '?', cluster))) %>%
  arrange(study, position)

### Analyse clusters and save results ###
plots <- make_cluster_plots(dms_wide, cols = !!cols, chem_env_cols = within_10_0_A:within_10_0_Y)
write_tsv(select(dms_wide, cluster, study, gene, position, wt), str_c('data/clusterings/', root_name, '.tsv'))
save_plotlist(plots, root = root_fig_dir, overwrite = 'all')
