#!/usr/bin/env Rscript
# Functions for importing studies in a standardised manner

source('src/config.R')

## general study data saving function
# dm_data = tibble with columns position, wt, mut, score, raw_score
# study_id = authour_year_gene style standard study id
# transform = string describing the transform applied
# fill = column to colour distribution plots by. If NULL no colouring is applied
standardise_study <- function(dm_data, study_id, transform = 'No Transform'){
  study_name = format_study(study_id)
  
  p_orig <- ggplot(dm_data, aes(x = raw_score, fill = class)) +
    guides(fill = guide_legend(title = 'Variant Class')) +
    scale_fill_manual(values = c(Missense='cornflowerblue', Nonsense='firebrick2', Synonymous='green2')) +
    geom_histogram(bins = 30) +
    labs(title = str_c('Original score distribution for ', study_name), x = 'Score', y = 'Count')

  p_trans <- ggplot(dm_data, aes(x = score, fill = class)) +
    guides(fill = guide_legend(title = 'Variant Class')) +
    scale_fill_manual(values = c(Missense='cornflowerblue', Nonsense='firebrick2', Synonymous='green2')) +
    geom_histogram(bins = 30) +
    labs(title = str_c('Transformed score distribution for ', study_name), x = 'Score', y = 'Count')
  
  # Write output
  if (!dir.exists(str_c('figures/0_data_properties/', study_id))){
    dir.create(str_c('figures/0_data_properties/', study_id))
  }
  ggsave(str_c('figures/0_data_properties/', study_id, '/original_distribution.pdf'), p_orig, units = 'cm', height = 12, width = 20)
  ggsave(str_c('figures/0_data_properties/', study_id, '/transformed_distribution.pdf'), p_trans, units = 'cm', height = 12, width = 20)
  
  write_tsv(select(dm_data, position, wt, mut, score), str_c('data/studies/', study_id, '/', study_id, '.tsv'))
}

## Function to determine variant class
get_variant_class <- function(wt, mut){
  if (!length(wt) == length(mut)){
    stop('wt and mut vectors must be the same length')
  }
  
  out <- rep('Missense', length(wt))
  out[wt == mut] <- 'Synonymous'
  out[mut == '*'] <- 'Nonsense'
  
  return(out)
}

## Scale and normalise VAMP-seq style
# data ranges from ~0 (NULL) -> 1 (wt) -) >1 beneficial
transform_vamp_seq <- function(x){
  # Transform
  y <- 1 + (x - 1) / -min(x - 1, na.rm = TRUE)
  y <- log2(y + min(y[y > 0], na.rm = TRUE))
  
  # Normalise
  y <- y / -min(y, na.rm = TRUE)
}