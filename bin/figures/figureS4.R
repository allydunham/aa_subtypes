#!/usr/bin/env Rscript
# Produce figure S4 (Study Confounding)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- read_tsv('data/combined_mutational_scans.tsv') %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))

umap2_breaks <- c(-2.5, 0, 2.5)

study_pretty <- sapply(unique(dms$study), format_study, max_width = 20)

p_confounding <- mutate(dms, study_pretty = study_pretty[study]) %>%
  ggplot(aes(x = umap1, y = umap2, colour = study_pretty)) +
  facet_wrap(~study_pretty, ncol = 4) +
  geom_point(data = dms, colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_y_continuous(breaks = umap2_breaks) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = FALSE)
ggsave('figures/4_figures/figureS4.pdf', p_confounding, width = 183, height = 270, units = 'mm')
ggsave('figures/4_figures/figureS4.png', p_confounding, width = 183, height = 270, units = 'mm')
