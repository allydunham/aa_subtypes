#!/usr/bin/env Rscript
# Assemble and Analyse the final set of subtypes, combining different values of deepSplit for different AAs
source('src/config.R')
source('src/subtype_characterisation.R')

### Setup ###
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

# Calculate Standard plots
plots <- plot_cluster_diagnostics(dms, clusters, cols = PC2:PC20)
plots <- c(plots, plot_cluster_characterisation(full_characterisation, selective_characterisation, clusters))

### Subtypes Heatmaps ###
# Outliers
outlier_profiles <- filter(dms, str_detect(cluster, str_c("^", CLUSTER_OUTLIER_RE, "$"))) %>%
  mutate(id = str_c(gene, position, sep = ' ')) %>%
  select(id, wt, A:Y) %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
  mutate(er = clamp(er, 2, -2),
         id = as.factor(id),
         mut = add_markdown(mut, AA_COLOURS))

plots$outlier_profiles <- (ggplot(outlier_profiles, aes(x = mut, y = id, fill = er)) +
                             lemon::facet_rep_grid(rows=vars(wt), space='free_y', scales = 'free', repeat.tick.labels = TRUE) +
                             geom_raster() +
                             scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction) +
                             labs(caption = str_wrap('Note: outliers (|ER| > 2) have been clamped, affecting a few positions near to 2 and two extreme values (|ER| > 4)', width = 60)) +
                             theme(axis.ticks = element_blank(),
                                   axis.text.x = element_markdown(),
                                   axis.title = element_blank(),
                                   strip.placement = 'outside',
                                   strip.text.y = element_text(angle = 0),
                                   panel.grid.major.y = element_blank())) %>%
  labeled_plot(units = 'cm', width = 15, height = 100)

# Main Subtypes
plot_cluster <- function(tbl, cluster, breaks){
  tbl <- mutate(tbl, id = str_c(gene, position, sep = ' ')) %>%
    select(id, A:Y) %>%
    pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
    mutate(mut = add_markdown(mut, AA_COLOURS))
  
  (ggplot(tbl, aes(x = mut, y = id, fill = er)) +
     geom_raster() +
     scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                          limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
     labs(title = cluster,
          caption = str_wrap('Note: outliers (|ER| > 1.5) have been clamped, mainly affecting a small number of extreme values (|ER| > 3)', width = 60)) +
     theme(axis.ticks = element_blank(),
           axis.text.x = element_markdown(),
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

### Frequency of permissive/not proline subtypes ###
not_proline_subtypes <- c('A3', 'D3', 'E2', 'G4', 'I3', 'K3', 'L6', 'M2', 'N2', 'Q2', 'R2', 'S2', 'T2', 'V5', 'Y4')
most_selective_subtypes <- group_by(full_characterisation$profiles, cluster) %>%
  summarise(mean_er = mean(er)) %>%
  mutate(aa = str_sub(cluster, end = 1)) %>%
  group_by(aa) %>%
  filter(!mean_er > min(mean_er))

freq_summary <- group_by(full_characterisation$summary, aa) %>%
  mutate(freq = n / sum(n)) %>%
  summarise(Permissive = freq[which(cluster == str_c(aa, 'P'))],
            `Not Proline` = max(freq[which(cluster %in% not_proline_subtypes)], 0),
            `Most Selective` = freq[which(cluster %in% most_selective_subtypes$cluster)],
            Other = 1 - (Permissive + `Not Proline` + `Most Selective`)) %>%
  pivot_longer(-aa, names_to = 'type', values_to = 'freq') %>%
  mutate(type = factor(type, levels = c('Permissive', 'Not Proline', 'Other', 'Most Selective'))) %>%
  left_join(select(most_selective_subtypes, aa, mean_er), by = 'aa') %>%
  mutate(aa = add_markdown(aa, AA_COLOURS))
er_limits <- pretty_break(most_selective_subtypes$mean_er, rough_n = 4, sym = 0, sig_figs = 2)

plots$subtype_freqs <- ggplot(freq_summary) +
  geom_col(aes(x = freq, y = aa, fill = type)) +
  geom_point(aes(x = -0.05, y = aa, colour = mean_er), shape = 15, size = 6) +
  scale_fill_brewer(type = 'qual', palette = 'Paired') + 
  scale_colour_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                       limits = er_limits$limits, breaks = er_limits$breaks, labels = er_limits$labels) +
  guides(fill = guide_legend(title = 'Subtype'),
         colour = guide_colourbar(title = str_wrap('Mean ER of Most Selective Subtype', width = 15))) +
  labs(x = 'Frequency') +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(),
        axis.ticks.y = element_blank())

### G1 / G3 Differences ###
g1_g3_terms <- filter(dms, cluster %in% c('G1', 'G3')) %>%
  select(cluster, solvation_polar, solvation_hydrophobic, van_der_waals, van_der_waals_clashes, entropy_sidechain,
         entropy_mainchain, torsional_clash, backbone_clash, phi, psi, starts_with('angstroms_to')) %>%
  mutate(Nearest_AA_Dist = apply(select(., angstroms_to_A:angstroms_to_Y), 1, function(x){min(x[x>3.9])})) %>% # shorter than this are neighbours in the chain
  select(-starts_with('angstroms')) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'value') %>%
  mutate(term = str_replace_all(term, '_', ' ') %>% str_to_title() %>% str_replace('Aa', 'AA'))

plots$g1_vs_g3 <- (ggplot(g1_g3_terms, aes(x = value, y = ..scaled.., colour = cluster)) +
                     facet_wrap(~term, scales = 'free_x', labeller = label_wrap_gen(20), strip.position = 'bottom', nrow = 3) +
                     stat_density(geom = 'line', position = 'identity') +
                     labs(y = 'Density') +
                     scale_colour_brewer(type = 'qual', palette = 'Set1') + 
                     guides(colour = guide_legend(title = '', override.aes = list(shape = 15))) +
                     theme(strip.placement = 'outside',
                           axis.title.x = element_blank())) %>%
  labeled_plot(width = 20, height = 20, units = 'cm')

### Large / Small hydrophobic Differences ###
large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
hydro_groups <- structure(rep(c('Large', 'Small', 'Aromatic'), each=3), names = c(large_hydrophobics, small_hydrophobics, aromatics))

hydrophobic_terms <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  select(cluster, solvation_polar, solvation_hydrophobic, van_der_waals, van_der_waals_clashes, entropy_sidechain,
         entropy_mainchain, torsional_clash, backbone_clash, phi, psi, side_chain_abs, starts_with('angstroms_to')) %>%
  mutate(Nearest_AA_Dist = apply(select(., angstroms_to_A:angstroms_to_Y), 1, function(x){min(x[x>3.9])})) %>% # shorter than this are neighbours in the chain
  select(-starts_with('angstroms')) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'value') %>%
  mutate(term = str_replace_all(term, '_', ' ') %>% str_to_title() %>% str_replace('Aa', 'AA'),
         group = hydro_groups[cluster])

plots$hydrophobic_size_difference <- (ggplot(hydrophobic_terms, aes(x = value, y = ..scaled.., colour = group)) +
                                        facet_wrap(~term, scales = 'free_x', labeller = label_wrap_gen(20), strip.position = 'bottom', nrow = 3) +
                                        stat_density(geom = 'line', position = 'identity') +
                                        labs(y = 'Density') +
                                        scale_colour_brewer(type = 'qual', palette = 'Dark2') + 
                                        guides(colour = guide_legend(title = '', override.aes = list(shape = 15))) +
                                        theme(strip.placement = 'outside',
                                              axis.title.x = element_blank())) %>%
  labeled_plot(width = 20, height = 20, units = 'cm')

### Save Results ###
write_tsv(select(dms, cluster, study, gene, position, wt), 'data/subtypes/final_subtypes.tsv')
saveRDS(clusters, file = 'data/subtypes/final_subtypes.rds')
save_plotlist(plots, root = 'figures/2_subtypes/final_subtypes/', overwrite = 'all')
