#!/usr/bin/env Rscript
# Assemble and Analyse the final set of subtypes, combining different values of deepSplit for different AAs
source('src/config.R')
source('src/subtype_characterisation.R')

#### Setup Data ####
config <- read_yaml('meta/final_subtypes.yaml')

# Import base data
dms <- read_tsv('data/combined_mutational_scans.tsv')
sections <- map(read_yaml('meta/structures.yaml'), extract2, 'sections')
pdb_pos <- dms_pdb_positions(dms, sections)
dms <- mutate(dms, pdb_position = pdb_pos$position, pdb_chain = pdb_pos$chain)

# Import clusterings
clusters_ds0 <- readRDS('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.rds')
clusters_ds1 <-readRDS('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.rds')
clusters <- matrix(c(clusters_ds0, clusters_ds1), ncol = 2) %>% set_rownames(names(clusters_ds0))

cluster_tbl <- bind_rows(`0`=read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.tsv'),
                         `1`=read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.tsv'),
                         .id='deepSplit') %>%
  mutate(deepSplit = as.integer(deepSplit))

# Assemble choosen clusters
aa_deep_split_inds <- unlist(config$deepSplit)[rownames(clusters)]

clusters <- clusters[cbind(1:20, aa_deep_split_inds + 1)] %>% set_names(names(config$deepSplit))
cluster_tbl <- left_join(tibble(wt=names(aa_deep_split_inds), deepSplit=aa_deep_split_inds),
                         cluster_tbl, by = c('wt', 'deepSplit'))

dms <- full_join(cluster_tbl, dms, by = c('study', 'gene', 'position', 'wt'))

# Calculate characterisation
full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

outlier_clusters <- filter(full_characterisation$profiles, str_detect(cluster, CLUSTER_PERMISSIVE_RE) | str_detect(cluster, CLUSTER_OUTLIER_RE)) %>%
  pull(cluster) %>%
  unique()
selective_characterisation <- full_cluster_characterisation(filter(dms, !cluster %in% outlier_clusters))
n_clusters_selective <- nrow(full_characterisation$summary)
########

#### Make standard analysis plots ####
plots <- plot_cluster_diagnostics(dms, clusters, cols = PC2:PC20)

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
plots$multi_position_subtype_consistency <- labeled_plot(plot_cluster_multiple_experiment_consistency(full_characterisation), units='cm', height = 20, width=20)

if (n_clusters_selective > 1){
  plots$er_correlation_selective <- labeled_plot(plot_cluster_profile_correlation(selective_characterisation), units='cm', width=0.25*n_clusters_selective + 2, height=0.25*n_clusters_selective + 2)
  plots$er_cosine_selective <- labeled_plot(plot_cluster_profile_cosine_sim(selective_characterisation), units='cm', width=0.25*n_clusters_selective + 2, height=0.25*n_clusters_selective + 2)
} else {
  plots$er_correlation_selective <- ggplot()
  plots$er_cosine_selective <- ggplot()
}

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

if (n_clusters_selective > 1){
  profiles_selective <- cluster_mean_profiles(filter(dms, !cluster %in% outlier_clusters)) 
  plots$overall_dend_selective <- labeled_plot(plot_profile_dendogram(profiles_selective, A:Y, distance_method = 'cosine'), width=40, height=20)
} else {
  plots$overall_dend_selective <- ggplot()
}

grouped_profiles <- mutate(profiles, aa = str_sub(cluster, end = 1)) %>%
  group_by(aa)
plots$aa_dends <- group_map(grouped_profiles, ~plot_profile_dendogram(., A:Y, distance_method = 'cosine')) %>% 
  set_names(group_keys(grouped_profiles)$aa)
########

#### Additional analysis for final subtypes ####
## Plot Subtypes Heatmaps
# Outliers
outlier_profiles <- filter(dms, str_detect(cluster, str_c("^", CLUSTER_OUTLIER_RE, "$"))) %>%
  mutate(id = str_c(gene, position, sep = ' ')) %>%
  select(id, wt, A:Y) %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
  mutate(er = clamp(er, 2, -2),
         id = as.factor(id))

plots$outlier_profiles <- (ggplot(outlier_profiles, aes(x = mut, y = id, fill = er)) +
                             lemon::facet_rep_grid(rows=vars(wt), space='free_y', scales = 'free', repeat.tick.labels = TRUE) +
                             geom_raster() +
                             scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction) +
                             labs(caption = str_wrap('Note: outliers (|ER| > 2) have been clamped, affecting a few positions near to 2 and two extreme values (|ER| > 4)', width = 60)) +
                             theme(axis.ticks = element_blank(),
                                   axis.text.x = element_text(colour = AA_COLOURS[sort(unique(outlier_profiles$mut))]),
                                   axis.title = element_blank(),
                                   strip.placement = 'outside',
                                   strip.text.y = element_text(angle = 0),
                                   panel.grid.major.y = element_blank())) %>%
  labeled_plot(units = 'cm', width = 15, height = 100)

# Main Subtypes
plot_cluster <- function(tbl, cluster, breaks){
  tbl <- mutate(tbl, id = str_c(gene, position, sep = ' ')) %>%
    select(id, A:Y) %>%
    pivot_longer(A:Y, names_to = 'mut', values_to = 'er')
  
  (ggplot(tbl, aes(x = mut, y = id, fill = er)) +
     geom_raster() +
     scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                          limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
     labs(title = cluster,
          caption = str_wrap('Note: outliers (|ER| > 1.5) have been clamped, mainly affecting a small number of extreme values (|ER| > 3)', width = 60)) +
     theme(axis.ticks = element_blank(),
           axis.text.x = element_text(colour = AA_COLOURS[sort(unique(outlier_profiles$mut))]),
           axis.title = element_blank(),
           strip.placement = 'outside',
           strip.text.y = element_text(angle = 0),
           panel.grid.major.y = element_blank(),
           plot.title = element_text(hjust = 0.5))) %>%
    labeled_plot(units = 'cm', width = 20, height = 30)
}

cluster_positions <- select(dms, cluster, gene, position, wt, A:Y) %>%
  filter(!str_detect(cluster, str_c("^", CLUSTER_OUTLIER_RE, "$"))) %>%
  mutate_at(vars(A:Y), ~clamp(., 1.5, -1.5)) %>%
  group_by(cluster)

breaks <- pivot_longer(cluster_positions, A:Y, names_to = 'mut', values_to = 'er') %>% 
  pull(er) %>% 
  pretty_break(rough_n = 3, sym = 0)

plots$cluster_heatmaps <- group_map(cluster_positions, ~plot_cluster(., cluster = .y, breaks = breaks)) %>%
  set_names(group_keys(cluster_positions)$cluster)

### Save Results ###
write_tsv(select(dms, cluster, study, gene, position, wt), 'data/subtypes/final_subtypes.tsv')
saveRDS(clusters, file = 'data/subtypes/final_subtypes.rds')
save_plotlist(plots, root = 'figures/2_subtypes/final_subtypes/', overwrite = 'all')
