#!/ussr/bin/env Rscript
# Perform a clustering of AA subtypes, based on a YAML config
source('src/config.R')
source('src/clustering.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Make and analyse AA subtypes', name = 'AA Subtypes')
parser <- add_argument(parser, arg = 'yaml', help = 'YAML file parameterising clustering')
parser <- add_argument(parser, arg = '--data', help = 'Directory to save cluster assignments', default = 'data/clustering')
parser <- add_argument(parser, arg = '--figures', help = 'Root directory to save figures', default = 'figures/2_clustering')
args <- parse_args(parser)

conf <- read_yaml(args$yaml)
root_name <- str_split(basename(args$yaml), '\\.', simplify = TRUE)[1]

cols <- eval(parse(text = str_c('quo(', conf$columns,')')))

dms_wide <- read_tsv('data/combined_mutational_scans.tsv')

### Setup function to pass arguments to cluster function ###
# need to use ... construct 
make_cluster_func <- function(func){
  return(function(tbl, ...){func(tbl, !!cols, ...)})
}

cluster_funcs <- list(kmeans=make_kmeans_clusters, hclust=make_hclust_clusters, hclust_dynamic=make_dynamic_hclust_clusters,
                      dbscan=make_dbscan_clusters, hdbscan=make_hdbscan_clusters)

if (conf$method %in% names(cluster_funcs)){
  cluster_func <- make_cluster_func(cluster_funcs[[conf$method]])
} else {
  stop(str_c('Unrecognised method "', conf$method, '", please use one of: ', str_c(names(cluster_funcs), collapse = ', ')))
}

### Make Clusters ###
clusters <- group_by(dms_wide, wt) %>%
  group_map(~do.call(cluster_func, c(list(tbl=.), conf$args)), keep = TRUE) %>%
  set_names(sapply(., function(x){first(x$tbl$wt)}))

dms_wide <- map_dfr(clusters, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, cluster)) %>%
  arrange(study, position)

### Analyse clusters and save results ###
plots <- make_cluster_plots(dms_wide, cols = !!cols, chem_env_cols = within_10_0_A:within_10_0_Y, clusters = clusters)

write_tsv(select(dms_wide, cluster, study, gene, position, wt), str_c(args$data, '/', root_name, '.tsv'))

root_fig_dir <- str_c(args$figures, '/', root_name)
dir.create(root_fig_dir, recursive = TRUE)
save_plotlist(plots, root = root_fig_dir, overwrite = 'all')
