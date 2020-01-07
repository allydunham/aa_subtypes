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

# Split numeric vector into limits/break points/labels inteligently for display on legends etc.
# Using rough_n is approximate
# Sym allows you to make the range symetrical around a given value (usually 0)
pretty_break <- function(x, step=NULL, rough_n=NULL, sig_figs=4, sym=NULL){
  if (is.null(step) & is.null(rough_n)){
    rough_n <- 3
  }
  
  limits <- range(x)
  if (!is.null(sym)){
    m <- max(abs(limits - sym))
    limits <- m * c(-1, 1) + sym
  }
  
  if (is.null(step) & !is.null(rough_n)){
    step = signif((limits[2] - limits[1])/rough_n, 1)
  }
  
  limits_trunc <- trunc(limits/step) * step
  breaks <- seq(limits_trunc[1], limits_trunc[2], step)
  if (abs(breaks[1]) < abs(limits[1])){
    breaks <- c(limits[1], breaks)
  }
  if (abs(breaks[length(breaks)]) < abs(limits[2])){
    breaks <- c(breaks, limits[2])
  }
  
  labels <- str_remove_all(signif(breaks, sig_figs), "\\.0*$")
  
  if (abs(breaks[1] - breaks[2]) < 0.5 * step){
    labels[2] <- ''
  }
  if (abs(breaks[length(breaks)] - breaks[length(breaks) - 1]) < 0.5 * step){
    labels[length(breaks) - 1] <- ''
  }
  
  return(list(limits=limits, breaks=breaks, labels=labels))
}

# Clamp a value between two limts
clamp <- function(x, upper=Inf, lower=-Inf){
  x[x > upper] <- upper
  x[x < lower] <- lower
  return(x)
}