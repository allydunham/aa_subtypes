#!/usr/bin/env Rscript
# Produce figure S5 (Proportion of variance)
source('src/config.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

pca <- tibble_pca(dms, A:Y)

pca_summary <- tibble(pc=0:20, sd=c(0, pca$sdev)) %>%
  mutate(prop_var = sd^2/sum(sd^2),
         cum_var = cumsum(prop_var))

figure <- ggplot(pca_summary, aes(x = pc)) +
  geom_col(aes(y = prop_var, fill = 'Explained\nVariance')) + 
  geom_line(aes(y = cum_var, colour = 'Cumulative\nExplained\nVariance')) +
  scale_fill_manual(values = c(`Explained\nVariance`='cornflowerblue'), name = '') +
  scale_colour_manual(values = c(`Cumulative\nExplained\nVariance`='red'), name = '') +
  labs(x = 'Principal Component', y = 'Proportion of Variance')
ggsave('figures/4_figures/figureS5.pdf', figure, width = 183, height = 100, units = 'mm')
ggsave('figures/4_figures/figureS5.png', figure, width = 183, height = 100, units = 'mm')
