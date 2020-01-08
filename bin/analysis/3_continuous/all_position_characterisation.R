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
