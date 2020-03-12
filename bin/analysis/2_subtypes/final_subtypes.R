#!/usr/bin/env Rscript
# Assemble and Analyse the final set of subtypes, combining different values of deepSplit for different AAs
source('src/config.R')
source('src/subtype_characterisation.R')
plots <- list()

### Import Data ###
config <- read_yaml('meta/final_subtypes.yaml')

# Import base data
dms <- read_tsv('data/combined_mutational_scans.tsv')
sections <- map(read_yaml('meta/structures.yaml'), extract2, 'sections')
pdb_pos <- dms_pdb_positions(dms, sections)
dms <- mutate(dms, pdb_position = pdb_pos$position, pdb_chain = pdb_pos$chain)

# Import clusterings
clusters_ds0 <- readRDS('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.rds')
clusters_ds1 <- readRDS('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.rds')

cluster_tbl_ds0 <- read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.tsv')
cluster_tbl_ds1 <- read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.tsv')

# Assemble choosen clusters


# Calculate characterisation
full_characterisation <- full_cluster_characterisation(select(dms, cluster = cluster_ds0, everything()))

### Analyse Outliers ###
outlier_profiles <- filter(dms, str_detect(cluster_ds0, str_c("^", CLUSTER_OUTLIER_RE, "$"))) %>%
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

### Plot clusters ###
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

cluster_positions <- select(dms, cluster=cluster_ds0, gene, position, wt, A:Y) %>%
  filter(!str_detect(cluster, str_c("^", CLUSTER_OUTLIER_RE, "$"))) %>%
  mutate_at(vars(A:Y), ~clamp(., 1.5, -1.5)) %>%
  group_by(cluster)

breaks <- pivot_longer(cluster_positions, A:Y, names_to = 'mut', values_to = 'er') %>% 
  pull(er) %>% 
  pretty_break(rough_n = 3, sym = 0)

plots$cluster_heatmaps <- group_map(cluster_positions, ~plot_cluster(., cluster = .y, breaks = breaks)) %>%
  set_names(group_keys(cluster_positions)$cluster)

### Save figures ###
save_plotlist(plots, root = 'figures/2_subtypes/final_subtypes/', overwrite = 'all')
