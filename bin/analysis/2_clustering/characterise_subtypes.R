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
plots <- list()

### Calculate Profiles for all clusters ###
full_characterisation <- full_cluster_characterisation(dms)

get_aa_plot <- function(x, global_scale=TRUE){
  clusters <- full_characterisation$summary$cluster[str_starts(full_characterisation$summary$cluster, x)]
  plot_full_characterisation(clusters, full_characterisation, exclude_outliers = TRUE, global_scale = global_scale)
}

### Plot Per AA characterisations ###
# aa_profiles and aa_profiles_relative refer to the scales being shared among all AAs or specific to each
plots$aa_profiles <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=TRUE)
plots$aa_profiles_relative <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=FALSE)

### Plot all-subtype summaries ###
plots$er_vs_surface_accessibility <- (filter(full_characterisation$summary, !str_ends(cluster, '0')) %>%
                                        ggplot(aes(x = mean_er, y = mean_sa, label = cluster, colour = aa)) +
                                        facet_wrap(~aa) + 
                                        geom_text() +
                                        geom_smooth(method = 'lm', se = FALSE, fullrange = TRUE) + 
                                        scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                                        labs(x = 'Mean Norm. ER'), y = 'Mean Surface Accessibility') %>%
  labeled_plot(width = 15, height = 10, units = 'cm')

plots$er_vs_size <- (filter(full_characterisation$summary, !str_ends(cluster, '0')) %>%
                       ggplot(aes(x = mean_er, y = n, label = cluster, colour = aa)) +
                       facet_wrap(~aa) + 
                       geom_text() +
                       scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                       labs(x = 'Mean Norm. ER'), y = 'Count') %>%
  labeled_plot(width = 15, height = 10, units = 'cm')

### Save Plots ###
root <- str_c(args$figures, '/', basename(file_path_sans_ext(args$subtypes)))
save_plotlist(plots, root, verbose = 3)


