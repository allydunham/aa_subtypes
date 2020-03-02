#!/usr/bin/env Rscript
# Analyse the final set of chosen subtypes
source('src/config.R')
source('src/subtype_characterisation.R')
plots <- list()

### Import Data ###
sections <- map(read_yaml('meta/structures.yaml'), extract2, 'sections')
dms <- left_join(rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0.tsv'), cluster_ds0 = cluster),
                 rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1.tsv'), cluster_ds1 = cluster),
                 by = c("study", "gene", "position", "wt")) %>%
  left_join(read_tsv('data/combined_mutational_scans.tsv'), by = c("study", "gene", "position", "wt")) %>%
  select(cluster_ds0, cluster_ds1, everything())

pdb_pos <- dms_pdb_positions(dms, sections)
dms <- mutate(dms, pdb_position = pdb_pos$position, pdb_chain = pdb_pos$chain)

full_characterisation <- full_cluster_characterisation(select(dms, cluster = cluster_ds0, everything()))

### Analyse Outliers ###
outlier_profiles <- filter(dms, str_detect(cluster_ds0, '^[A-Z]0$')) %>%
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
                                   panel.grid.major.y = element_blank(),
                                   plot.title = element_text(hjust = 0))) %>%
  labeled_plot(units = 'cm', width = 15, height = 100)

### Save figures ###
save_plotlist(plots, root = 'figures/2_subtypes/final_subtypes/', overwrite = 'all')
