#!/usr/bin/env Rscript
# Produce figure 6 (Large/Small Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

er_limits <- c(min(full_characterisation$profiles$er), -min(full_characterisation$profiles$er))

plot_profiles <- function(clusters, legend = FALSE){
  filter(full_characterisation$profiles, cluster %in% clusters) %>%
    mutate(mut = add_markdown(mut, AA_COLOURS),
           cluster = add_markdown(cluster, cluster_number_colourmap(cluster))) %>%
    ggplot(aes(x = mut, y = cluster, fill = er)) +
    geom_raster(show.legend = legend) +
    coord_fixed() +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
    guides(fill = guide_colourbar(title = 'Normalised ER')) + 
    theme(axis.text.x = element_markdown(),
          axis.text.y = element_markdown(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major.y = element_blank())
}

large_hydrophobics <- c('I2', 'L2', 'M1')
small_hydrophobics <- c('A1', 'G2', 'P3')
aromatics <- c('F2', 'W1', 'Y1')
hydro_groups <- structure(rep(c('Large Aliphatic', 'Small Aliphatic', 'Aromatic'), each=3), names = c(large_hydrophobics, small_hydrophobics, aromatics))

hydro_profiles <- filter(full_characterisation$profiles, cluster %in% names(hydro_groups)) %>%
  mutate(group = hydro_groups[cluster]) %>%
  group_by(group, mut) %>%
  summarise(er = mean(er)) %>%
  ungroup() %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         group = as.factor(group))

hydro_group_labs <- add_markdown(levels(hydro_profiles$group), c(`Large Aliphatic`='#D95F02', `Small Aliphatic`='#7570B3', Aromatic='#1B9E77')) %>%
  str_c('**', ., '**')
hydro_subtype_labs <- map(levels(hydro_profiles$group),
                          ~names(hydro_groups[hydro_groups == .]) %>% 
                            add_markdown(., cluster_number_colourmap(.)) %>%
                            str_c(collapse = ', ')) %>%
  set_names(hydro_group_labs) %>% unlist()

p_hydro_profiles <- ggplot(hydro_profiles, aes(x = mut, y = as.integer(group), fill = er)) +
  geom_raster(show.legend = FALSE) +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(hydro_group_labs), labels = hydro_group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(hydro_subtype_labs), labels = hydro_subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER')) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.text.y.left = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_hydro_sa <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  mutate(group = hydro_groups[cluster]) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = group)) +
  geom_density(show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  labs(x = 'Surface Accessibility (All Atom Abs.)', y = 'Scaled Density') +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

p_hydro_foldx <- filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
  select(cluster, entropy_mainchain, entropy_sidechain, solvation_polar, solvation_hydrophobic,
         van_der_waals, van_der_waals_clashes, torsional_clash, backbone_clash) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'ddg') %>%
  mutate(group = hydro_groups[cluster],
         term_pretty = FOLDX_TERMS[term]) %>%
  ggplot(aes(x = group, y = ddg, fill = group)) +
  facet_wrap(~term_pretty, scales = 'free_x', nrow = 2) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill = FALSE) +
  labs(x = '', y = expression('Mean FoldX Prediction (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

plot_hydro_foldx_boxplot <- function(term){
  term <- enquo(term)
  
  filter(dms, cluster %in% c(large_hydrophobics, small_hydrophobics, aromatics)) %>%
    mutate(group = hydro_groups[cluster]) %>%
    ggplot(aes(x = group, y = !!term, fill = group)) +
    geom_boxplot(show.legend = FALSE) +
    geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
    scale_fill_brewer(type = 'qual', palette = 'Dark2') +
    labs(x = '', y = parse(text = str_c(FOLDX_TERMS_PLOTMATH[quo_name(term)],"~Delta*Delta*'G (kj mol'^-1*')'"))) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.placement = 'outside',
          legend.position = 'bottom')
}

p_hydro_foldx_entropy_sidechain <- plot_hydro_foldx_boxplot(entropy_sidechain)
p_hydro_foldx_solvation_hydrophobic <- plot_hydro_foldx_boxplot(solvation_hydrophobic)
p_hydro_foldx_van_der_waals <- plot_hydro_foldx_boxplot(van_der_waals)
p_hydro_foldx_van_der_waals_clashes <- plot_hydro_foldx_boxplot(van_der_waals_clashes)
p_hydro_foldx_backbone_clash <- plot_hydro_foldx_boxplot(backbone_clash)

### Assemble figure ###
size <- theme(text = element_text(size = 8))
er_legend <- plot_profiles('A1', legend = TRUE) %>% get_legend() %>% as_ggplot() + size

figure6 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size + 
  labs(tag = 'D') +
  annotation_custom(ggplotGrob(p_hydro_profiles + size), xmin = 0, xmax = 1, ymin = 0.38, ymax = 0.62) +
  annotation_custom(ggplotGrob(p_hydro_foldx_entropy_sidechain + size), xmin = 0, xmax = 0.25, ymin = 0.65, ymax = 1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ras_aliphatic_entropy.png'), interpolate = TRUE,
                    xmin = 0.275, xmax = 0.475, ymin = 0.75, ymax = 0.95) +
  annotate('text', x = 0.25, y = 0.65, label = str_wrap('Small residues are more entropically favourable', 40)) +
  annotation_custom(ggplotGrob(p_hydro_foldx_van_der_waals_clashes + size), xmin = 0.5, xmax = 0.75, ymin = 0.65, ymax = 1) +
  annotation_raster(readPNG('figures/4_figures/position_examples/adrb2_ala_small_hydro.png'), interpolate = TRUE,
                    xmin = 0.775, xmax = 0.975, ymin = 0.75, ymax = 0.95) +
  annotate('text', x = 0.75, y = 0.65, label = str_wrap('Small hydrophobic positions are cramped, and larger substitutions lead to clashes', 40)) +
  annotation_custom(ggplotGrob(p_hydro_foldx_van_der_waals + size), xmin = 0, xmax = 0.25, ymin = 0, ymax = 0.35) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_phe_pi.png'), interpolate = TRUE,
                    xmin = 0.275, xmax = 0.475, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.25, y = 0.05, label = str_wrap('Larger residues, particularly aromatics, create stronger Van der Waals forces', 40)) +
  annotation_custom(ggplotGrob(p_hydro_foldx_solvation_hydrophobic + size), xmin = 0.5, xmax = 0.75, ymin = 0, ymax = 0.35) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ras_met_buried.png'), interpolate = TRUE,
                    xmin = 0.775, xmax = 0.975, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.75, y = 0.05, label = str_wrap('Larger residues have a bigger solvation energy contribution', 40))

ggsave('figures/4_figures/figure6.pdf', figure6, width = 200, height = 200, units = 'mm')
ggsave('figures/4_figures/figure6.png', figure6, width = 200, height = 200, units = 'mm')
