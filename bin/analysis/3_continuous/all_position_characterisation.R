#!/usr/bin/env Rscript
# Characterise based on continuous ER gradiant
source('src/config.R')
source('src/continuous.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

### Correlations between various factors and Mean ER ###
mean_er_cors <- select(dms, wt, mean_score, backbone_hbond:energy_ionisation, all_atom_abs, within_10_0_A:within_10_0_Y, ss_g:ss_t) %>%
  group_by(wt) %>%
  group_modify(~as_tibble(cor(select(., mean_score), select(., -mean_score), use = "pairwise.complete.obs"))) %>%
  ungroup() %>%
  pivot_longer(-wt, names_to = 'term', values_to = 'cor') %>%
  mutate(type = get_factor_type(term), term = str_remove(term, 'ss_|within_10_0_')) %>%
  add_factor_order(wt, term, cor)

term_labels <- c(structure(LETTERS, names=LETTERS), all_atom_abs="'All Atom Abs.'", FOLDX_TERMS_PLOTMATH, DSSP_CLASSES_PLOTMATH)
p_er_cors <- ggplot(mean_er_cors, aes(x = wt, y = term, fill = cor)) +
  facet_grid(rows = vars(type), scales = 'free_y', space = 'free_y', drop = TRUE) +
  geom_raster() +
  scale_fill_gradient2(guide = guide_colourbar(title = 'Pearson\nCorrelation')) +
  scale_y_discrete(labels = sapply(term_labels, function(x){parse(text=x)})) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = AA_COLOURS[levels(mean_er_cors$wt)]))
ggsave('figures/3_continuous/mean_er_correlations.pdf', p_er_cors, units = 'cm', width = 25, height = 30)

### Per AA PCA ###
aa_pcas <- group_by(dms, wt)
aa_pcas <- group_map(aa_pcas, ~tibble_pca(., A:Y)) %>%
  set_names(group_keys(aa_pcas)$wt)

## PC Loadings
aa_pca_rotations <- map(aa_pcas, ~as_tibble(t(.$rotation), rownames = 'pc')) %>%
  bind_rows(.id = 'wt') %>%
  pivot_longer(A:Y, names_to = 'mut', values_to = 'rotation') %>%
  mutate(pc = factor(pc, levels = str_c('PC', 20:1)))

plot_aa_pca_loadings <- function(x, ...){
  breaks <- pretty_break(x$rotation, rough_n = 5, sym = 0)
  (ggplot(x, aes(x = mut, y = pc, fill = rotation)) +
      geom_raster() +
      scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                           limits = breaks$limits, breaks=breaks$breaks, labels=breaks$labels) +
      coord_fixed() +
      guides(fill=guide_colourbar(title = 'Rotation')) +
      theme(axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(colour = AA_COLOURS[sort(unique(x$mut))]))) %>%
    labeled_plot(units = 'cm', width = 10, height = 10)
}

p_aa_pca <- group_by(aa_pca_rotations, wt)
p_aa_pca <- group_map(p_aa_pca, plot_aa_pca_loadings) %>%
  set_names(group_keys(p_aa_pca)$wt)
save_plotlist(p_aa_pca, root = 'figures/3_continuous/per_aa_pca_loadings')  

## PC levels
dms_aa_pca <- map(unique(dms$wt), ~bind_cols(filter(dms, wt == .), as_tibble(aa_pcas[[.]]$x) %>% rename_all(~str_c('aa', .)))) %>%
  bind_rows() %>%
  arrange(study, position)

