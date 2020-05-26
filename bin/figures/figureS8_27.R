#!/usr/bin/env Rscript
# Produce figure S8-27 (Subtype Characterisation)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)

figures <- group_by(dms, wt) %>%
  group_map(~plot_full_characterisation(unique(.$cluster), full_characterisation, exclude_outliers = TRUE, global_scale = FALSE)) %>%
  map(extract2, 'overall') %>%
  set_names(str_c('figureS', 8:27))

save_plotlist(figures, 'figures/4_figures/', default_format = 'pdf')
save_plotlist(figures, 'figures/4_figures/', default_format = 'png')

# Save single PDF version
pdf('figures/4_figures/figureS8_27.pdf', onefile = TRUE, width = 11.7, height = 8.3)
for (name in names(figures)){
  p <- ggplot() +
    geom_blank() +
    lims(x=c(0,1), y=c(0,1)) +
    labs(caption = str_c('Figure ', str_sub(name, start = -2))) +
    annotation_custom(figures[[name]], xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    theme(plot.caption = element_text(hjust = 0.5),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(8, 1, 1, 1), 'mm'))
  print(p)
}
dev.off()
