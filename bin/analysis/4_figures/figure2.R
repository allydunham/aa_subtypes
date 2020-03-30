#!/usr/bin/env Rscript
# Produce figure 2 (Mutational Landscape)
source('src/config.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

### Panel 1 - Mean Score ###
p_mean_score <- ggplot(dms, aes(x = PC1, y = mean_score)) +
  geom_point() +
  labs(x = 'PC1', y = 'Mean Normalised ER')

### Panel 2 - Surface Accessibility ###
p_surface_accessibility <- drop_na(dms, side_chain_rel) %>%
  ggplot(aes(x = umap1, y = umap2, colour = all_atom_abs)) +
  geom_point() +
  scale_colour_gradientn(colours = c('#1a2a6c', '#b21f1f', '#fdbb2d'),
                         values = c(0, 0.2, 1)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = str_wrap('Surface Accessibility', 10)))

### Panel 3 - AA hydrophobicity ### 
p_hydrophobicity <- drop_na(dms, hydrophobicity) %>%
  ggplot(aes(x = umap1, y = umap2, colour = hydrophobicity)) +
  geom_point() +
  scale_colour_gradientn(colours = c('#4575b4', '#e0f3f8', '#fee090', '#fc8d59', '#d73027'),
                         values = rescale01(c(-0.4, 0, 0.4, 0.8, 1.2))) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'Hydrophobicity'))

### Panel 4 - Similar Domains ###
p_similar_domains <- ggplot(tibble(x = 1, y = 1), aes(x=x, y=y)) +
  geom_point(colour = 'white') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank()) +
  annotate('text', x = 1, y = 1, label = 'Similar domains\nmapped to UMAP')

### Assemble Figure ###
size <- theme(text = element_text(size = 10))
p1 <- p_mean_score + labs(tag = 'A') + size
p2 <- p_surface_accessibility + labs(tag = 'B') + size
p3 <- p_hydrophobicity + labs(tag = 'C') + size
p4 <- p_similar_domains + labs(tag = 'D') + size

figure2 <- multi_panel_figure(width = 300, height = 200, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)
ggsave('figures/4_figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')
ggsave('figures/4_figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')

