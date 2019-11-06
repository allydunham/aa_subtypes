#!/usr/bin/env Rscript
# Test whether studies on the same gene relate well to each other

source('src/config.R')

# Import a study
import_study <- function(d){
  study <- str_split(d, '/')[[1]]
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
  coord_equal()

#### HSP90 ####
########

#### TEM1 ####
########

#### UBI ####
########