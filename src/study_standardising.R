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
  
  write_tsv(select(dm_data, position, wt, mut, score, raw_score, class), str_c('data/studies/', study_id, '/', study_id, '.tsv'))
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

# Calculate E-score equivalent to Enrich1 
# Currently no pseudocount etc. (simple implementation without error checking)
e_score <- function(sel, bkg){
  bkg[bkg == 0] <- NA
  
  freq_sel <- sel/sum(sel, na.rm = TRUE)
  freq_bkg <- bkg/sum(bkg, na.rm = TRUE)
  
  return(freq_sel/freq_bkg)
}

#### Functions for Melnikov et al. 2014 (APH(3')-II) ####
# Read aa count tables from melnikov et al. 2014
read_melnikov_table <- function(fi){
  tbl <- read_tsv(str_c('data/studies/melnikov_2014_aph3ii/raw/', fi), skip = 1, col_names = FALSE, col_types = cols(.default = col_character())) %>%
    t() %>%
    as_tibble(rownames = NULL, .name_repair='minimal') %>%
    set_colnames(.[1,]) %>%
    filter(!Position == 'Position') %>%
    rename(position = Position,
           wt = `Wild-type`) %>%
    mutate_at(vars(-wt), as.numeric)
  return(tbl)
}

# Wrapper to pass correct background and selection counts to fitness function, based on format of Melnikov 2014 data
# Expects sel to be a data.frame with cols for position, wt and all mut's in one selection/drug/library category
# these are given as exp_name in the SX_DRUG_LX format of melnikov
melnikov_fitness <- function(sel, exp_name, bkg){
  # Extract meta info on experiment
  meta <- as.list(strsplit(exp_name, '_')[[1]])
  names(meta) <- c('selection_round', 'drug', 'library')
  
  # Select correct background reads for library
  bkg <- bkg[[str_c('Bkg', str_sub(meta$library, -1))]]
  
  # Format bkg and sel as matrices
  ref_aas <- bkg$wt
  gene_length <- length(ref_aas)
  sel <- as.matrix(select(sel, -position, -wt))
  bkg <- as.matrix(select(bkg, -position, -wt))
  
  # Apply simple pseudocount of minimum non zero
  pseudo <- min(sel[sel>0], na.rm = TRUE)
  sel <- sel + pseudo
  bkg <- bkg + pseudo
  
  # Calculate e-score per position row - this allows calculation of ER for each variant taking account of other positions
  # as that information is contained in the positional wt count
  e_scores <- t(sapply(1:nrow(sel), function(x){e_score(sel[x,], bkg[x,])}))
  
  # Not properly possible to tell how fully WT sequence fairs as the WT AA measures include lots of mutants too
  # So cannot normalise to WT, however the per position method does leave most WT residues at ~1 already so scale stands
  fitness <- log2(e_scores + min(e_scores[e_scores > 0], na.rm = TRUE)) %>%
    as_tibble(.name_repair = 'unique') %>%
    mutate(position = 1:gene_length,
           wt = ref_aas) %>%
    gather(key = 'mut', value = 'score', -wt, -position)
  return(fitness)
}
########