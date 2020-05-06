#!/usr/bin/env Rscript
# Produce figure S2 (Gene Repeats)
source('src/config.R')
source('src/study_standardising.R')

### BRCA1 ###
brca1_studys <- sapply(c('data/studies/findlay_2018_brca1', 'data/studies/starita_2015_brca1'), import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score) %>% 
  drop_na()

p_brca1 <- ggplot(brca1_studys, aes(x=findlay_2018_brca1, y=starita_2015_brca1, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  labs(x = 'Findlay et al. 2018', y = 'Starita et al. 2015', title = 'BRCA1') +
  guides(colour = FALSE) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

### HSP90 ###
# No overlap between hietpas/jiang and mishra studies
hsp90_studys <- sapply(c('data/studies/hietpas_2011_hsp90', 'data/studies/jiang_2013_hsp90'),
                       import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_hsp90 <- ggplot(hsp90_studys, aes(x=hietpas_2011_hsp90, y=jiang_2013_hsp90, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  labs(x = 'Hietpas et al. 2011', y = 'Jiang et al. 2013', title = 'HSP90') +
  guides(colour = FALSE) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

### TEM1 ###
## Both same lab
tem1_studys <- sapply(c('data/studies/firnberg_2014_tem1', 'data/studies/steinberg_2016_tem1'),
                      import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_tem1 <- ggplot(tem1_studys, aes(x=firnberg_2014_tem1, y=steinberg_2016_tem1, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  labs(x = 'Firnberg et al. 2014', y = 'Steinberg & Ostermeier 2016', title = 'TEM1') +
  guides(colour = FALSE) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

# Ubi
ubi_studys <- sapply(c('data/studies/roscoe_2013_ubi', 'data/studies/roscoe_2014_ubi'),
                     import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_ubi <- ggplot(ubi_studys, aes(x=roscoe_2013_ubi, y=roscoe_2014_ubi, colour=class)) +
  geom_point(shape=20) +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  scale_x_continuous(expand = expansion(0.01)) +
  scale_y_continuous(expand = expansion(0.01)) +
  guides(colour = guide_legend(title = '', direction = 'horizontal')) +
  labs(x = 'Roscoe et al. 2013', y = 'Roscoe & Bolon 2014', title = 'UBI') +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())

p_legend <- get_legend(p_ubi) %>% as_ggplot()
p_ubi <- p_ubi + guides(colour = FALSE)

# Assemble
figure <- multi_panel_figure(width = c(89, 89), height = c(5, 89, 89), unit = 'mm',
                             row_spacing = 0, column_spacing = 0, panel_label_type = 'none') %>%
  fill_panel(p_legend, row = 1, column = 1:2) %>%
  fill_panel(p_brca1 + labs(tag = 'A'), row = 2, column = 1) %>%
  fill_panel(p_hsp90 + labs(tag = 'B'), row = 2, column = 2) %>%
  fill_panel(p_tem1 + labs(tag = 'C'), row = 3, column = 1) %>%
  fill_panel(p_ubi + labs(tag = 'D'), row = 3, column = 2)
ggsave('figures/4_figures/figureS2.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS2.png', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

