#!/usr/bin/env Rscript
# Utility functions for Subtypes project

# Function to generate a pretty study name from the id using YAML config data
format_study <- function(x, max_width=60, study_dir='data/studies', mark_filtered=FALSE){
  yaml <- read_yaml(str_c(study_dir, '/', x, '/', x, '.yaml'))
  
  study <- str_c(yaml$authour, ' ', yaml$year, ' (', yaml$gene, ')') %>%
    str_wrap(width = max_width)
  
  if (mark_filtered & yaml$qc$filter){
    study <- str_c(study, '*')
  }
  
  return(study)
}

# Consistent filename from gene
gene_to_filename <- function(x){
  return(str_replace_all(str_to_lower(x), ' ', ''))
}

# Make pretty p value categories for plots
pretty_p_values <- function(p, breaks = c(0.001, 0.01, 0.05)){
  p_out <- rep(str_c('> ', breaks[length(breaks)]), length(p))
  for (cutoff in sort(breaks, decreasing = TRUE)){
    p_out[p < cutoff] <- str_c('< ', cutoff)
  }
  p_out <- factor(p_out, levels = c(str_c('> ', breaks[length(breaks)]),
                                    sapply(sort(breaks, decreasing = TRUE), function(x){str_c('< ', x)})))
  return(p_out)
}

# Transform standardised data table to wide format
make_dms_wide <- function(dms){
  foldx_averages <- select(dms, study, position, wt, total_energy:entropy_complex) %>%
    select(-sloop_entropy, -mloop_entropy, -entropy_complex, -water_bridge) %>% # Drop terms that are unused in our structures
    drop_na(total_energy) %>%
    group_by(study, position, wt) %>%
    summarise_all(mean, na.rm=TRUE)
  
  position_constants <- select(dms, study, position, wt, phi:hydrophobicity) %>%
    distinct()
  
  dms_wide <- filter(dms, mut %in% Biostrings::AA_STANDARD) %>%
    select(study, gene, position, wt, mut, imputed_score, log10_sift) %>%
    pivot_wider(names_from = mut, values_from = c(imputed_score, log10_sift)) %>%
    rename_at(vars(starts_with('imputed_score_')), ~str_sub(., start=-1))

  mutate(dms_wide,
         mean_score = rowMeans(select(dms_wide, A:Y)),
         mean_sift = rowMeans(select(dms_wide, log10_sift_A:log10_sift_Y))) %>%
    left_join(foldx_averages, by = c('study', 'position', 'wt')) %>%
    left_join(position_constants, by = c('study', 'position', 'wt'))
  
}

