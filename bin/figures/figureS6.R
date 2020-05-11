#!/usr/bin/env Rscript
# Produce figure S6 (Frequency of rarer subtypes)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

sub_freqs <- group_by(full_characterisation$summary, aa) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  filter(!str_detect(cluster, CLUSTER_PERMISSIVE_RE), !str_detect(cluster, CLUSTER_OUTLIER_RE)) %>% 
  mutate(num = as.integer(str_sub(cluster, 2)))

model <- lm(log(freq) ~ num + 0, data = sub_freqs)
coefs <- summary(model)$coefficients['num',]
x <- c(seq(1, 8, 0.1), seq(8, 1, -0.1))
y <- c(exp((coefs['Estimate']+coefs['Std. Error'])*seq(1, 8, 0.1)),
       exp((coefs['Estimate']-coefs['Std. Error'])*seq(8, 1, -0.1)))

figure <- ggplot(sub_freqs, aes(x = num, y = freq)) + 
  geom_point(aes(colour = aa)) +
  annotate('polygon', x = x, y = y, fill = 'lightgrey', alpha = 0.5, colour = NA) +
  stat_function(fun = ~exp(coefs['Estimate']*.)) +
  scale_colour_manual(values = AA_COLOURS) +
  scale_x_continuous(breaks = 1:8) +
   guides(colour = guide_legend(title = '', direction = 'horizontal', nrow = 2, byrow = TRUE)) +
   labs(x = 'Subtype', y = 'Frequency') +
   theme(legend.position = 'bottom')

ggsave('figures/4_figures/figureS6.pdf', figure, width = 187, height = 100, units = 'mm')
ggsave('figures/4_figures/figureS6.png', figure, width = 187, height = 100, units = 'mm')

