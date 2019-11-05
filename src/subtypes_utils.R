#!/usr/bin/env Rscript
# Utility functions for Subtypes project

# Function to generate a pretty study name from the id using YAML config data
format_study <- function(x, max_width=60, study_dir='data/studies'){
  yaml <- read_yaml(str_c(study_dir, '/', x, '/', x, '.yaml'))
  
  str_c(yaml$authour, ' ', yaml$year, ' (', yaml$gene, ')') %>%
    str_wrap(width = max_width) %>%
    return()
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
