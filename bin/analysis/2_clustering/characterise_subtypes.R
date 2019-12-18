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

get_aa_plot <- function(x, global_scale=TRUE){
  clusters <- full_characterisation$summary$cluster[str_starts(full_characterisation$summary$cluster, x)]
  plot_full_characterisation(clusters, full_characterisation, exclude_outliers = TRUE, global_scale = global_scale)
}

# global and relative refer to the scales being shared among all AAs or specific to each
plots_global <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=TRUE)
plots_relative <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=FALSE)

root_global <- str_c(args$figures, '/', basename(file_path_sans_ext(args$subtypes)), '/aa_profiles')
root_relative <- str_c(args$figures, '/', basename(file_path_sans_ext(args$subtypes)), '/aa_profiles_relative')

dir.create(root_global)
dir.create(root_relative)

for (aa in names(plots_global)){
  message(aa)
  save_multi_panel_figure(plots_global[[aa]]$overall, str_c(root_global, '/', aa, '.pdf'))
}

for (aa in names(plots_relative)){
  message(aa)
  save_multi_panel_figure(plots_relative[[aa]]$overall, str_c(root_relative, '/', aa, '.pdf'))
}
