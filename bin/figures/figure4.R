#!/usr/bin/env Rscript
# Produce figure 4 (Subtype Examples)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

### Panel 1 - C examples ###
p_cys <- blank_plot('C examples')

### Panel 2 - D/E examples ###
p_charge <- blank_plot('D/E (R/K) Examples')

### Panel 3 - G1/G3 examples ###
p_gly <- blank_plot('G1/G3 Example')

### Panel 4 - Large/Small hydrophobics ###
p_hydro <- blank_plot('Large/Small/Aromatic Example')

### Assemble figure ###
size <- theme(text = element_text(size = 8))
p1 <- p_cys + labs(tag = 'A') + size
p2 <- p_charge + labs(tag = 'B') + size
p3 <- p_gly + labs(tag = 'C') + size
p4 <- p_hydro + labs(tag = 'D') + size

figure4 <- multi_panel_figure(width = 200, height = 200, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)
ggsave('figures/4_figures/figure4.pdf', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')
ggsave('figures/4_figures/figure4.png', figure4, width = figure_width(figure4), height = figure_height(figure4), units = 'mm')

