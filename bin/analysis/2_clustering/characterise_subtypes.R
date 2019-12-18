#!/usr/bin/env Rscript
# Characterise generated clusters
source('src/config.R')
source('src/clustering.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Characterise AA subtypes', name = 'AA Subtype Characerisation')
parser <- add_argument(parser, arg = 'subtypes', help = 'TSV file assigning positions to clusters')
parser <- add_argument(parser, arg = '--dms', help = 'Path to DMS data', default = 'data/combined_mutational_scans.tsv')
parser <- add_argument(parser, arg = '--figures', help = 'Root directory to save figures', default = 'figures/2_clustering')
args <- parse_args(parser)

subtypes <- read_tsv(args$subtypes)
dms <- read_tsv(args$dms) %>%
  left_join(subtypes, ., by = c('study', 'gene', 'position', 'wt'))

full_characterisation <- full_cluster_characterisation(dms)

get_aa_plot <- function(x){
  clusters <- full_characterisation$summary$cluster[str_starts(full_characterisation$summary$cluster, x)]
  plot_full_characterisation(clusters, full_characterisation, exclude_outliers = TRUE)
}
plots <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE)

root <- str_c(args$figures, '/', basename(file_path_sans_ext(args$subtypes)),'/aa_profiles')
dir.create(root)

for (aa in names(plots)){
  message(aa)
  save_multi_panel_figure(plots[[aa]]$overall, str_c(root, '/', aa, '.pdf'))
}
