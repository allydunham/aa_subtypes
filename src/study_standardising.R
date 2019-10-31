#!/usr/bin/env Rscript
# Functions for importing studies in a standardised manner

source('src/config.R')

AA_THREE_2_ONE <- structure(names(Biostrings::AMINO_ACID_CODE), names = Biostrings::AMINO_ACID_CODE)
AA_THREE_2_ONE['Ter'] <- '*'

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

## Normalise Score
normalise_score <- function(x){
  x / -min(x, na.rm = TRUE)
}

## Scale and normalise VAMP-seq style
# data ranges from ~0 (NULL) -> 1 (wt) -) >1 beneficial
transform_vamp_seq <- function(x, normalise=TRUE){
  # Transform
  y <- 1 + (x - 1) / -min(x - 1, na.rm = TRUE)
  y <- log2(y + min(y[y > 0], na.rm = TRUE))
  
  # Normalise
  if (normalise){
    return(normalise_score(y))
  } else {
    return(y)
  }
  
}

# Calculate E-score equivalent to Enrich1 
# Currently no pseudocount etc. (simple implementation without error checking)
e_score <- function(sel, bkg){
  bkg[bkg == 0] <- NA
  
  freq_sel <- sel/sum(sel, na.rm = TRUE)
  freq_bkg <- bkg/sum(bkg, na.rm = TRUE)
  
  return(freq_sel/freq_bkg)
}

## Import MAVEDB study
read_mavedb <- function(path, score_col=NULL, score_transform=identity, position_offset = 0){
  score_col <- enquo(score_col)
  if (rlang::quo_is_null(score_col)){
    score_col <- quo(score)
  }
  
  read_csv(path, skip = 4) %>%
    tidyr::extract(hgvs_pro, into = c('wt', 'position', 'mut'), "p.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})", convert = TRUE) %>%
    mutate(wt = AA_THREE_2_ONE[wt], mut = AA_THREE_2_ONE[mut], position = position + position_offset) %>%
    rename(raw_score = !!score_col) %>%
    mutate(score = normalise_score(score_transform(raw_score)),
           class = get_variant_class(wt, mut)) %>%
    select(position, wt, mut, score, raw_score, class) %>%
    arrange(position, mut) %>%
    return()
}

# Untangle seqIDs of the form 1,2-A,D
process_split_seqid <- function(x){
  x <- str_split(x, '[-,]')[[1]]
  return(str_c(x[1:(length(x)/2)], x[(length(x)/2 + 1):length(x)], collapse = ','))
}

# Get muts from seq, expects each as a character vector
muts_from_seq <- function(mut_seq, wt_seq){
  if (all(mut_seq == wt_seq)){
    return(NA)
  }
  
  pos <- which(!mut_seq == wt_seq)
  return(str_c(wt_seq[pos], pos, mut_seq[pos], collapse = ','))
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

#### Functions for Kitzman et al. 2015 (GAL4) ####
read_kitzman_sheet <- function(path, sheet){
  tbl <- read_xlsx(path, skip = 1, na = 'ND', sheet = sheet) %>%
    rename(position = `Residue #`) %>%
    mutate(wt = apply(., 1, function(x, nam){nam[x == 'wt' & !is.na(x)]}, nam = names(.)),
           label = sheet) %>%
    gather(key = 'mut', value = 'log2_enrichment', -position, -wt, -label) %>%
    mutate(log2_enrichment = if_else(log2_enrichment == 'wt', '0', log2_enrichment)) %>% # set wt to 0 log2 enrichment ratio
    mutate(log2_enrichment = as.numeric(log2_enrichment))
  return(tbl)
}
########

#### Functions for Mishra et al. 2016 (HSP90) ####
read_mishra_sheet <- function(path, sheet){
  tbl <- read_xlsx(path, sheet = sheet, col_names = str_c('col', 1:13))
  
  # Check sheet type
  if (tbl[1,1] == 'Stop counts'){
    ## Process sheets with a single replicate
    nom <- tbl[7,] %>% unlist(., use.names = FALSE)
    tbl <- tbl[8:nrow(tbl),] %>%
      set_names(nom) %>%
      rename_at(vars(-position, -aa), list( ~ paste0('rep1_', .))) %>%
      mutate_at(vars(-aa), as.numeric) %>%
      mutate(avg = rep1_norm_ratiochange)
    
  } else {
    ## Process sheets with replicates
    # Get first row of sub-tables
    top_row <- which(tbl$col1 == 'position') + 1
    
    # Get bottom row of sub-tables
    bot_row <- c(which(apply(tbl, 1, function(x){all(is.na(x))})) - 1, dim(tbl)[1])
    bot_row <- sapply(top_row, function(x){bot_row[which(bot_row > x)[1]]})
    
    # Extract sub-table names
    rep_nom <- tbl[top_row[1] - 1,] %>% unlist(., use.names = FALSE)
    ave_nom <- tbl[top_row[length(top_row)] - 1,] %>% unlist(., use.names = FALSE)
    ave_nom <- ave_nom[!is.na(ave_nom)]
    
    # Extract Subtables and add names
    rep1 <- tbl[top_row[1]:bot_row[1],] %>% 
      set_names(rep_nom) %>%
      rename_at(vars(-position, -aa), list( ~ paste0('rep1_', .)))
    
    rep2 <- tbl[top_row[2]:bot_row[2],] %>%
      set_names(rep_nom) %>%
      rename_at(vars(-position, -aa), list( ~ paste0('rep2_', .)))
    
    ave <- tbl[top_row[3]:bot_row[3],] %>%
      select_if(colSums(!is.na(.)) > 0) %>%
      set_names(ave_nom) %>%
      select(-s1, -s2) %>% # also found in rep tbls
      rename(aa = `amino acid`)
    
    tbl <- full_join(ave, rep1,, by=c('position', 'aa')) %>%
      full_join(., rep2, by=c('position', 'aa')) %>%
      mutate_at(vars(-aa), as.numeric)
  }
  return(tbl)
}
########