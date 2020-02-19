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
plots <- list()

### Calculate Profiles for all clusters ###
full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

### Plot all cluster characterisation ###
plots$ramachandran_angles <- labeled_plot(plot_cluster_ramachandran_angles(full_characterisation), units='cm', width=20, height=20)
plots$sizes <- labeled_plot(plot_cluster_sizes(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$er_profiles <- labeled_plot(plot_cluster_profiles(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$er_correlation <- labeled_plot(plot_cluster_profile_correlation(full_characterisation), units='cm', width=0.25*n_clusters + 2, height=0.25*n_clusters + 2)
plots$er_distance <- labeled_plot(plot_cluster_profile_distances(full_characterisation), units='cm', width=0.25*n_clusters + 2, height=0.25*n_clusters + 2)
plots$er_cosine <- labeled_plot(plot_cluster_profile_cosine_sim(full_characterisation), units='cm', width=0.25*n_clusters + 2, height=0.25*n_clusters + 2)
plots$foldx <- labeled_plot(plot_cluster_foldx_profiles(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$chem_env <- labeled_plot(plot_cluster_chem_env_profiles(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$aa_distance <- labeled_plot(plot_cluster_aa_distances(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$ss_probability <- labeled_plot(plot_cluster_ss_profile(full_characterisation), units='cm', width=0.75*n_clusters + 2, height=20)
plots$profile_variance <- group_by(full_characterisation$tbl, wt)
plots$profile_variance <- group_map(plots$profile_variance, ~labeled_plot(plot_cluster_profile_variation(.), units='cm', height=20, width=30), keep = TRUE) %>%
  set_names(group_keys(plots$profile_variance)$wt)
plots$multi_position_subtype_consistency <- labeled_plot(plot_cluster_multiple_experiment_consistency(x), units='cm', height = 20, width=20)

### Plot Per AA characterisations ###
get_aa_plot <- function(x, global_scale=TRUE){
  clusters <- full_characterisation$summary$cluster[str_starts(full_characterisation$summary$cluster, x)]
  plot_full_characterisation(clusters, full_characterisation, exclude_outliers = TRUE, global_scale = global_scale)
}

# relative and global refer to the scales being shared among all AAs or specific to each
plots_global <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=TRUE)
plots_relative <- sapply(unique(str_sub(full_characterisation$summary$cluster, end = 1)), get_aa_plot, simplify = FALSE, global_scale=FALSE)

plots$aa_profiles <- map(plots_global, extract2, 'overall')
plots$aa_profiles_relative <- map(plots_relative, extract2, 'overall')

### Plot all-subtype summaries ###
plots$er_vs_surface_accessibility <- (filter(full_characterisation$summary, !str_detect(cluster, '^[A-Z]0$')) %>%
                                        ggplot(aes(x = mean_er, y = mean_sa, label = cluster, colour = aa)) +
                                        facet_wrap(~aa) + 
                                        geom_text() +
                                        geom_smooth(method = 'lm', se = FALSE, fullrange = TRUE) + 
                                        scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                                        labs(x = 'Mean Norm. ER', y = 'Mean Surface Accessibility')) %>%
  labeled_plot(width = 20, height = 15, units = 'cm')

plots$er_vs_size <- (filter(full_characterisation$summary, !str_detect(cluster, '^[A-Z]0$')) %>%
                       ggplot(aes(x = mean_er, y = n, label = cluster, colour = aa)) +
                       facet_wrap(~aa) + 
                       geom_text() +
                       scale_colour_manual(values = AA_COLOURS, guide = FALSE) + 
                       labs(x = 'Mean Norm. ER', y = 'Count')) %>%
  labeled_plot(width = 20, height = 15, units = 'cm')

plots$ss_probabilities <- group_by(full_characterisation$tbl, wt) %>%
  group_map(~labeled_plot(plot_cluster_ss_density(.), units = 'cm', height = 20, width = 20), keep = TRUE) %>%
  set_names(sort(unique(full_characterisation$tbl$wt)))

### Plot compressed dendrograms for hclust ###
if ('hclust' %in% names(clusters[[1]])){
  plots$minimal_dends <- plot_compressed_dendrograms(clusters, dms)
}

### Plot profile dendograms ###
profiles <- cluster_mean_profiles(dms) 
plots$overall_dend <- labeled_plot(plot_profile_dendogram(profiles, A:Y, distance_method = 'cosine'), width=40, height=20)

grouped_profiles <- mutate(profiles, aa = str_sub(cluster, end = 1)) %>%
  group_by(aa)
plots$aa_dends <- group_map(grouped_profiles, ~plot_profile_dendogram(., A:Y, distance_method = 'cosine')) %>% 
  set_names(group_keys(grouped_profiles)$aa)

### Save Plots ###
save_plotlist(plots, args$figures, verbose = 2)
