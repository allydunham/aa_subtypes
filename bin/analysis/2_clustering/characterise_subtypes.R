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
# relative and global refer to the scales being shared among all AAs or specific to each
plots_global <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=TRUE)
plots_relative <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=FALSE)

plots$aa_profiles <- map(plots_global, extract2, 'overall')
plots$aa_profiles_relative <- map(plots_relative, extract2, 'overall')

### Plot all-subtype summaries ###
plots$er_vs_surface_accessibility <- (filter(full_characterisation$summary, !str_ends(cluster, '0')) %>%
                                        ggplot(aes(x = mean_er, y = mean_sa, label = cluster, colour = aa)) +
                                        facet_wrap(~aa) + 
                                        geom_text() +
                                        geom_smooth(method = 'lm', se = FALSE, fullrange = TRUE) + 
                                        scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                                        labs(x = 'Mean Norm. ER', y = 'Mean Surface Accessibility')) %>%
  labeled_plot(width = 20, height = 15, units = 'cm')

plots$er_vs_size <- (filter(full_characterisation$summary, !str_ends(cluster, '0')) %>%
                       ggplot(aes(x = mean_er, y = n, label = cluster, colour = aa)) +
                       facet_wrap(~aa) + 
                       geom_text() +
                       scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                       labs(x = 'Mean Norm. ER', y = 'Count')) %>%
  labeled_plot(width = 20, height = 15, units = 'cm')

plots$ss_probabilities <- filter(full_characterisation$tbl, !str_ends(cluster, '0')) %>%
  group_by(wt) %>%
  group_map(~labeled_plot(plot_ss_density(.), units = 'cm', height = 20, width = 20), keep = TRUE) %>%
  set_names(sort(unique(full_characterisation$tbl$wt)))

### Save Plots ###
root <- str_c(args$figures, '/', basename(file_path_sans_ext(args$subtypes)))
save_plotlist(plots, root, verbose = 2)
