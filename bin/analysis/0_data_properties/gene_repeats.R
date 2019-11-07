#!/usr/bin/env Rscript
# Test whether studies on the same gene relate well to each other

source('src/config.R')

dir.create('figures/0_data_properties/gene_repeats')

# Import a study
import_study <- function(d){
  study <- str_split(str_replace(d, '/$', ''), '/')[[1]]
  study <- study[length(study)]
  yaml = read_yaml(str_c(d, '/', study, '.yaml'))
  
  tbl <- read_tsv(str_c(d, '/', study, '.tsv')) %>%
    mutate(study = yaml$study,
           gene = yaml$gene)
  return(tbl)
}

#### BRCA1 ####
brca1_studys <- sapply(c('data/studies/findlay_2018_brca1', 'data/studies/starita_2015_brca1'), import_study, simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_brca <- ggplot(brca1_studys, aes(x=findlay_2018_brca1, y=starita_2015_brca1, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Findlay et al. 2018', y = 'Starita et al. 2015', title = 'Scores in BRCA1 from two studies') +
  guides(colour = guide_legend(title = 'Variant Type'))
ggsave('figures/0_data_properties/gene_repeats/brca1.pdf', p_brca, units = 'cm', width = 15, height = 17)

#### HSP90 ####
hsp90_studys <- sapply(c('data/studies/hietpas_2011_hsp90', 'data/studies/jiang_2013_hsp90', 'data/studies/mishra_2016_hsp90'),
                       import_study, simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

# No overlap between hietpas/jiang and mishra studies

p_hsp90 <- ggplot(hsp90_studys, aes(x=hietpas_2011_hsp90, y=jiang_2013_hsp90, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Hietpas et al. 2011', y = 'Jiang et al. 2013', title = 'Scores in HSP90 from two studies',
       subtitle = 'Both from the Bolon Lab') +
  guides(colour = guide_legend(title = 'Variant Type'))
ggsave('figures/0_data_properties/gene_repeats/hsp90.pdf', p_hsp90, units = 'cm', width = 10, height = 10)
########

#### TEM1 ####
## Both same lab
tem1_studys <- sapply(c('data/studies/firnberg_2014_tem1', 'data/studies/steinberg_2016_tem1'),
                       import_study, simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_tem1 <- ggplot(tem1_studys, aes(x=firnberg_2014_tem1, y=steinberg_2016_tem1, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Firnberg et al. 2014', y = 'Steinberg & Ostermeier 2016', title = 'Scores in TEM1 from two studies',
       subtitle = 'Both from the Ostermeier Lab') +
  guides(colour = guide_legend(title = 'Variant Type'))
ggsave('figures/0_data_properties/gene_repeats/tem1.pdf', p_tem1, units = 'cm', width = 10, height = 10)
########

#### UBI ####
ubi_studys <- sapply(c('data/studies/roscoe_2013_ubi', 'data/studies/roscoe_2014_ubi'),
                      import_study, simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_ubi <- ggplot(ubi_studys, aes(x=roscoe_2013_ubi, y=roscoe_2014_ubi, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Roscoe et al. 2013', y = 'Roscoe & Bolon 2014', title = 'Scores in UBI from two studies',
       subtitle = 'Both from the Bolon Lab') +
  guides(colour = guide_legend(title = 'Variant Type'))
ggsave('figures/0_data_properties/gene_repeats/ubi.pdf', p_ubi, units = 'cm', width = 10, height = 10)
########
