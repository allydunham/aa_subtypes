#!/usr/bin/env Rscript
# Produce figure 4 (Subtype Examples)
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

### Example 1 - C examples ###
p_cys_profiles <- plot_profiles(c('C1', 'C2'))

p_cys_surface_acc <- filter(dms, cluster %in% c('C1', 'C2')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density', show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_cys_disulphide <- filter(dms, cluster %in% c('C1', 'C2'), !is.na(disulfide)) %>%
  ggplot(aes(x = cluster, fill = cluster)) +
  geom_bar(show.legend = FALSE) +
  labs(y = 'Disulphide Bonds', x = '') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  theme(legend.title = element_blank(),
        panel.grid.major.y = element_blank())

### Example 2 - D/E examples ###
charge_groups <- c(D1='Negative', D2='Polar', D3='Not Proline', E1='Negative', E2='Not Proline', E3='Polar',
                   K1='Positive', K2='Not Negative', K3='Not Proline', R1='Positive', R2='Not Proline', R3='Not Negative')

charge_profiles <- filter(full_characterisation$profiles, cluster %in% names(charge_groups)) %>%
  mutate(group = charge_groups[cluster]) %>%
  group_by(group, mut) %>%
  summarise(er = mean(er)) %>%
  ungroup() %>%
  mutate(mut = add_markdown(mut, AA_COLOURS),
         group = as.factor(group))

charge_group_labs <- levels(charge_profiles$group)
charge_subtype_labs <- c(Negative="<span style = 'color:#e41a1c'>D1</span>, <span style = 'color:#e41a1c'>E1</span>",
                  Positive="<span style = 'color:#e41a1c'>K1</span>, <span style = 'color:#e41a1c'>R1</span>",
                  Polar="<span style = 'color:#377eb8'>D2</span>, <span style = 'color:#4daf4a'>E3</span>",
                  `Not Proline`="<span style = 'color:#4daf4a'>D3</span>, <span style = 'color:#377eb8'>E2</span>, <span style = 'color:#4daf4a'>K3</span>, <span style = 'color:#377eb8'>R2</span>",
                  `Not Negative`="<span style = 'color:#377eb8'>K2</span>, <span style = 'color:#4daf4a'>R3</span>")[charge_group_labs]

p_charge_profiles <- ggplot(charge_profiles, aes(x = mut, y = as.integer(group), fill = er)) +
  geom_raster(show.legend = FALSE) +
  coord_fixed() +
  scale_y_continuous(breaks = 1:length(charge_group_labs), labels = charge_group_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(charge_subtype_labs), labels = charge_subtype_labs)) +
  scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction, limits = er_limits) +
  guides(fill = guide_colourbar(title = 'Normalised ER')) + 
  theme(axis.text.x = element_markdown(),
        axis.text.y.right = element_markdown(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank())

p_charge_foldx <- filter(dms, cluster %in% names(charge_groups)) %>%
  mutate(group = charge_groups[cluster]) %>%
  ggplot(aes(x = group, y = electrostatics)) +
  geom_boxplot(show.legend = FALSE, fill = '#377eb8', outlier.shape = 20) +
  geom_hline(yintercept = 0, linetype = 'dotted', colour = 'black') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Subtype')) +
  labs(x = '', y = expression('Mean Substitution Electrostatic'~Delta*Delta*'G (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

### Example 3 - G1/G3 examples ###
p_gly_profiles <- plot_profiles(c('G1', 'G2', 'G3'))

p_gly_surface_acc <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density') +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

p_gly_foldx <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  select(cluster, entropy_sidechain, phi, psi, solvation_hydrophobic, solvation_polar, van_der_waals, van_der_waals_clashes) %>%
  pivot_longer(-cluster, names_to = 'term', values_to = 'ddg') %>%
  mutate(term_pretty = c(FOLDX_TERMS, phi='Phi', psi='Psi')[term]) %>%
  ggplot(aes(x = cluster, y = ddg, fill = cluster)) +
  facet_wrap(~term_pretty, scales = 'free_x', nrow = 2) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  guides(fill = FALSE) +
  labs(x = '', y = expression('Mean FoldX Prediction (kj mol'^-1*')')) +
  theme(panel.grid.major.y = element_blank())

p_gly_rama <- filter(dms, cluster %in% c('G1', 'G2', 'G3')) %>%
  ggplot(aes(x = phi, y = psi, colour = cluster)) +
  facet_wrap(~cluster) +
  geom_point() +
  coord_equal()

### Example 4 - Tyrosine ###
p_tyr_profiles <- plot_profiles(c('Y1', 'Y2', 'Y3', 'Y4', 'YP'))

p_tyr_surface_acc <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  ggplot(aes(x = all_atom_abs, y = ..scaled.., colour = cluster)) +
  geom_line(stat = 'density', show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Surface Accessibility (All Atom Abs)', y = 'Scaled Density') +
  theme(legend.title = element_blank())

tyr_dssp <- filter(dms, cluster %in% c('Y1', 'Y2', 'Y3', 'Y4', 'YP')) %>%
  select(cluster, study, position, starts_with('ss_')) %>%
  pivot_longer(starts_with('ss_'), names_to = 'term', names_prefix = 'ss_', values_to = 'p') %>%
  group_by(cluster, study, position) %>%
  summarise(ss = DSSP_CLASSES_PLOTMATH[term[which.max(p)]]) %>%
  ungroup() %>%
  count(cluster, ss) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

p_tyr_dssp <- ggplot(tyr_dssp, aes(x = cluster, y = prop, fill = ss)) +
  geom_col() +
  scale_fill_brewer(labels = function(x){parse(text=x)}, type = 'qual', palette = 'Dark2') +
  coord_flip() +
  labs(y = 'Proportion', x = '') +
  guides(fill = guide_legend(title = '', reverse = TRUE)) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')

p_tyr_haem <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/cbs_tyr_haem.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p_tyr_pi <- ggplot() +
  geom_blank() +
  annotation_custom(readPNG('figures/4_figures/position_examples/tp53_tyr_pi.png') %>% rasterGrob(interpolate = TRUE),
                    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

### Example 5 - Large/Small hydrophobics ###
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

# Panel 1
p1 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size +
  labs(tag = 'A') +
  annotation_custom(ggplotGrob(p_cys_profiles + size), xmin = 0.15, xmax = 0.8, ymin = 0.4, ymax = 0.6) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_disulphide.png'), interpolate = TRUE,
                    xmin = 0.025, xmax = 0.2, ymin = 0.7, ymax = 0.875) +
  annotation_custom(ggplotGrob(p_cys_disulphide + size), xmin = 0.2, xmax = 0.475, ymin = 0.7, ymax = 0.875) +
  annotate('text', x = 0.25, y = 0.65, label = str_wrap('Most disulphide bonds are C1 positions', 25)) + 
  annotation_custom(ggplotGrob(p_cys_surface_acc + size), xmin = 0.525, xmax = 0.975, ymin = 0.7, ymax = 0.975) +
  annotate('text', x = 0.75, y = 0.65, label = str_wrap('C2 positions are more buried', 25)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/gal4_cys_zinc.png'), interpolate = TRUE,
                    xmin = 0.1, xmax = 0.4, ymin = 0.1, ymax = 0.4) +
  annotate('text', x = 0.25, y = 0.05, label = str_wrap('Other active roles are generally C1 positions', 25)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_cys_aromatic.png'), interpolate = TRUE,
                    xmin = 0.525, xmax = 0.73, ymin = 0.15, ymax = 0.355) +
  annotation_raster(readPNG('figures/4_figures/position_examples/np_cys_aromatic.png'), interpolate = TRUE,
                    xmin = 0.77, xmax = 0.975, ymin = 0.15, ymax = 0.355) +
  annotate('text', x = 0.75, y = 0.05, label = str_wrap('Both C1 & C2 positions have aromatic interactions', 30))

# Panel 2
p2 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size + 
  labs(tag = 'B') +
  annotation_custom(ggplotGrob(p_charge_profiles + size), xmin = 0, xmax = 1, ymin = 0.35, ymax = 0.68) +
  annotate('text', x = 0.5, y = 0.975, label = str_wrap('Charged subtype positions have the strongest electrostatic interactions', 80)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_asp_ionic.png'), interpolate = TRUE,
                    xmin = 0.025, xmax = 0.225, ymin = 0.725, ymax = 0.925) +
  annotate('text', x = 0.125, y = 0.7, label = str_wrap('Ionic interactions', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_ligand.png'), interpolate = TRUE,
                    xmin = 0.275, xmax = 0.475, ymin = 0.725, ymax = 0.925) +
  annotate('text', x = 0.375, y = 0.7, label = str_wrap('Ligand binding', 30)) + 
  annotation_custom(ggplotGrob(p_charge_foldx + size), xmin = 0.5, xmax = 1, ymin = 0.65, ymax = 0.925) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tem1_asp_sa.png'), interpolate = TRUE,
                    xmin = 0, xmax = 0.25, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.125, y = 0.3, label = str_wrap('Many polar positions are surface accessible', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_arg_not_proline.png'), interpolate = TRUE,
                    xmin = 0.375, xmax = 0.625, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.5, y = 0.3, label = str_wrap('Not proline positions are often in secondary structures or tight turns', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_arg_not_neg.png'), interpolate = TRUE,
                    xmin = 0.75, xmax = 1, ymin = 0, ymax = 0.25) +
  annotate('text', x = 0.875, y = 0.3, label = str_wrap('Positions near negative residues or ligands (e.g. DNA backbone) need to not be negative them selves', 30))

# Panel 3
p3 <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted')) +
  size + 
  labs(tag = 'C') +
  annotation_custom(ggplotGrob(p_tyr_profiles + size), xmin = 0, xmax = 1, ymin = 0.38, ymax = 0.62) + 
  annotation_custom(ggplotGrob(p_tyr_surface_acc + size), xmin = 0.025, xmax = 0.475, ymin = 0.75, ymax = 1) +
  annotate('text', x = 0.25, y = 0.7, label = str_wrap('Permissive positions are surface accessible', 30)) +
  annotation_custom(ggplotGrob(p_tyr_dssp + size), xmin = 0.525, xmax = 0.975, ymin = 0.75, ymax = 1) + 
  annotation_raster(readPNG('figures/4_figures/position_examples/cbs_tyr_haem.png'), interpolate = TRUE,
                    xmin = 0.025, xmax = 0.225, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.125, y = 0.125, label = str_wrap('Y1 ligand binding', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/tp53_tyr_pi.png'), interpolate = TRUE,
                    xmin = 0.275, xmax = 0.475, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.375, y = 0.125, label = str_wrap('Pi orbital interactions', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/aph3ii_tyr_buried.png'), interpolate = TRUE,
                    xmin = 0.525, xmax = 0.725, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.625, y = 0.125, label = str_wrap('Y2 positions are buried', 30)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/pab1_tyr_not_proline.png'), interpolate = TRUE,
                    xmin = 0.775, xmax = 0.975, ymin = 0.15, ymax = 0.35) +
  annotate('text', x = 0.875, y = 0.125, label = str_wrap('Y4 maintaining backbone orientation', 30))
  
# Panel 4
p4 <- ggplot() +
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

# Final Figure
figure4 <- multi_panel_figure(width = 440, height = 400, columns = 11, rows = 10,
                              panel_label_type = 'none', row_spacing = 0.1, column_spacing = 0.1) %>%
  fill_panel(er_legend, row = 5:6, column = 11) %>%
  fill_panel(p1, row = 1:5, column = 1:5) %>%
  fill_panel(p2, row = 1:5, column = 6:10) %>%
  fill_panel(p3, row = 6:10, column = 1:5) %>%
  fill_panel(p4, row = 6:10, column = 6:10)

ggsave('figures/4_figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
