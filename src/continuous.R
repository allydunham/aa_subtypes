#!/usr/bin/env Rscript
# Functions for characterising Deep Mutational Scanning positions based on continuous ER gradiant

# Process the various factors
get_factor_type <- function(x){
  out <- rep(NA, length(x))
  out[x == 'all_atom_abs'] <- 'SA'
  out[x %in% names(FOLDX_TERMS)] <- 'FoldX'
  out[str_starts(x, 'ss_')] <- 'DSSP'
  out[str_starts(x, 'within_10_0')] <- 'Chem. Env.'
  
  return(out)
}

pretty_factors <- function(x){
  out <- rep(NA, length(x))
  typ <- get_factor_type(x)
  
  out[typ == 'SA'] <- 'All Atom Abs.'
  out[typ == 'FoldX'] <- FOLDX_TERMS[x[typ == 'FoldX']]
  out[typ == 'Secondary Structure'] <- DSSP_CLASSES_STR[str_sub(x[typ == 'Secondary Structure'], start = 4)]
  out[typ == 'Chemical Environment'] <- str_sub(x[typ == 'Chemical Environment'], start = -1)
    
  return(out)
}
