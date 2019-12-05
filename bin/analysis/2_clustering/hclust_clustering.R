#!/ussr/bin/env Rscript
# Perform an Hclust clustering of AA subtypes

source('src/config.R')
source('src/clustering.R')

### Parse Args and perform setup ###
library(argparser)
parser <- arg_parser(description = 'Make and analyse AA subtypes using kmeans clustering', name = 'Kmeans AA Clustering')
parser <- add_argument(parser, arg = '--height', help = 'Height to cut dendogram', default = NA, type = 'numeric')
parser <- add_argument(parser, arg = '--number', help = 'Number of clusters to select', default = NA, type = 'numeric')
parser <- add_argument(parser, arg = '--min_size', help = 'Minimum cluster size to consider', default = 5)
parser <- add_argument(parser, arg = '--mode', help = 'Cluster using "profile" or "pca"', default = 'profile')
parser <- add_argument(parser, arg = '--distance', help = 'Distance metric to use', default = 'manhattan')
args <- parse_args(parser)

if (!xor(is.na(args$number), is.na(args$height))){
  stop('Exactly one of --height and --number must be used')
}

if (!args$mode %in% names(CLUSTER_COLS)){
  stop(str_c('--mode must be one of ', str_c('"', names(CLUSTER_COLS), '"', collapse = ', ')))
}

if (!args$distance %in% c('euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski')){
  stop('--distance must be one of those supported by the R dist() function')
}

if (!is.na(args$height)){
  root_name <- str_c('hclust', args$mode, 'height', args$height, 'min', args$min_size, 'distance', args$distance, sep = '_')
} else if (!is.na(args$number)){
  root_name <- str_c('hclust', args$mode, 'number', args$number, 'min', args$min_size, 'distance', args$distance, sep = '_')
}

root_fig_dir <- str_c('figures/2_clustering/', root_name)
dir.create(root_fig_dir, recursive = TRUE)
dir.create('data/clusterings')

dms_wide <- read_tsv('data/combined_mutational_scans.tsv')

### Create Clusters ###
cols <- CLUSTER_COLS[[args$mode]]
hclust_cluster <- group_by(dms_wide, wt) %>%
  group_map(~make_hclust_clusters(., !!cols, h = args$height, k = args$number, min_size = args$min_size, dist_method=args$distance), keep = TRUE)

dms_wide <- map_dfr(hclust_cluster, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, cluster)) %>%
  arrange(study, position)

### Analyse clusters and save results ###
plots <- make_cluster_plots(dms_wide, cols = !!cols, chem_env_cols = within_10_0_A:within_10_0_Y)
write_tsv(select(dms_wide, cluster, study, gene, position, wt), str_c('data/clusterings/', root_name, '.tsv'))
save_plotlist(plots, root = root_fig_dir, overwrite = 'all')
