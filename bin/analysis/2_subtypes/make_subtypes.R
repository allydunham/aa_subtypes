#!/ussr/bin/env Rscript
# Perform a clustering of AA subtypes, based on a YAML config
source('src/config.R')
source('src/subtype_clustering.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Make and analyse AA subtypes', name = 'AA Subtypes')
parser <- add_argument(parser, arg = 'yaml', help = 'YAML file parameterising clustering')
parser <- add_argument(parser, arg = '--figures', help = 'Directory to save figures', default = NA)
parser <- add_argument(parser, arg = '--dms', help = 'DMS data path', default = 'data/combined_mutational_scans.tsv')
parser <- add_argument(parser, arg = '--out', help = 'Root path to output classifications', default = NA)
args <- parse_args(parser)

conf <- read_yaml(args$yaml)

if (is.na(args$out)){
  args$out <- str_split(basename(args$yaml), '\\.', simplify = TRUE)[1]
}

cols <- eval(parse(text = str_c('quo(', conf$columns,')')))

dms <- read_tsv(args$dms)

### Setup function to pass arguments to cluster function ###
# need to use ... construct 
make_cluster_func <- function(func){
  return(function(tbl, ...){func(tbl, !!cols, ...)})
}

cluster_funcs <- list(kmeans=make_kmeans_clusters, kmeans_multi=make_multi_kmeans_clusters,
                      hclust=make_hclust_clusters, hclust_dynamic=make_dynamic_hclust_clusters,
                      dbscan=make_dbscan_clusters, hdbscan=make_hdbscan_clusters, gmm=make_gmm_clusters, pam=make_pam_clusters)

if (conf$method %in% names(cluster_funcs)){
  cluster_func <- make_cluster_func(cluster_funcs[[conf$method]])
} else {
  stop(str_c('Unrecognised method "', conf$method, '", please use one of: ', str_c(names(cluster_funcs), collapse = ', ')))
}

### Make Clusters ###
clusters <- group_by(dms, wt) %>%
  group_map(~do.call(cluster_func, c(list(tbl=.), conf$args)), keep = TRUE) %>%
  set_names(sapply(., function(x){first(x$tbl$wt)}))

### Save Clusters ###
dms <- map_dfr(clusters, .f = ~ .$tbl) %>%
  mutate(cluster = str_c(wt, cluster)) %>%
  arrange(study, position)
write_tsv(select(dms, cluster, study, gene, position, wt), str_c(args$out, '.tsv'))
saveRDS(clusters, file = str_c(args$out, '.rds'))


### Save diagnostic plots ###
if (!is.na(args$figures)){
  plots <- plot_cluster_diagnostics(clusters, cols = !!cols)
  save_plotlist(plots, root = args$figures, overwrite = 'all')
}
