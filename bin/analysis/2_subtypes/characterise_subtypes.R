#!/usr/bin/env Rscript
# Characterise generated clusters
source('src/config.R')
source('src/subtype_characterisation.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Characterise AA subtypes', name = 'AA Subtype Characerisation')
parser <- add_argument(parser, arg = 'subtypes', help = 'Root filename assigning positions to subtypes. Should have a tsv (subtype per position) and rds (list of cluster objects) file.')
parser <- add_argument(parser, arg = '--dms', help = 'Path to DMS data', default = 'data/combined_mutational_scans.tsv')
parser <- add_argument(parser, arg = '--figures', help = 'Directory to save figures', default = '.')
args <- parse_args(parser)

subtypes <- read_tsv(str_c(args$subtypes, '.tsv'))
clusters <- readRDS(str_c(args$subtypes, '.rds'))

dms <- read_tsv(args$dms) %>%
  left_join(subtypes, ., by = c('study', 'gene', 'position', 'wt'))

### Calculate Profiles for all clusters ###
full_characterisation <- full_cluster_characterisation(dms)

# Make profiles with permissive/outlier clusters excluded
outlier_clusters <- filter(full_characterisation$profiles, str_detect(cluster, CLUSTER_PERMISSIVE_RE) | str_detect(cluster, CLUSTER_OUTLIER_RE)) %>%
  pull(cluster) %>%
  unique()
selective_characterisation <- full_cluster_characterisation(filter(dms, !cluster %in% outlier_clusters))
n_clusters_selective <- nrow(full_characterisation$summary)

### Plot all cluster characterisation ###
plots <- plot_cluster_characterisation(full_characterisation, selective_characterisation, clusters)

### Save Plots ###
save_plotlist(plots, args$figures, verbose = 2)
