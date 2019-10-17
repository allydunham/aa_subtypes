#!/usr/bin/env Rscript 
# Load libraries and configurations for all mutations analyses, to keep consistency within project

# Load all libraries (keeps namespace consistent)
library(Biostrings)
library(tidyverse)
library(rlang)
library(magrittr)
#library(broom)
library(readxl)
library(ggpubr)
library(GGally)
#library(scales)
#library(Rtsne)
library(yaml)

# Source custom functions
# source('src/util/misc_util.R')
# source('src/util/autoplot.R')
# source('src/data_processing/deep_mut_data.R')

#### Project config ####
# Categories of secondary structure
SS_REDUCED_HASH <- c(C='None', S='Turn', H='Helix', T='Turn', E='Strand', G='Helix', B='Strand', I='Helix')

# Categories of amino acid
# From Sigma-Aldrich website (https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
AA_REDUCED_CLASSES <- list(Aliphatic=c('A', 'I', 'L', 'M', 'V'),
                           Aromatic=c('F', 'W', 'Y'),
                           Polar=c('N', 'C', 'Q', 'S', 'T'),
                           Basic=c('R', 'H', 'K'),
                           Acidic=c('D', 'E'),
                           Glycine='G',
                           Proline='P')
AA_REDUCED_HASH <- structure(rep(names(AA_REDUCED_CLASSES), times=sapply(AA_REDUCED_CLASSES, length)),
                             names=unlist(AA_REDUCED_CLASSES))

## AA colours, based roughly on groups
AA_COLOURS <- c(A='red', I='salmon', L='firebrick', M='orange', V='tomato',
                F='gold', W='yellow3', Y='khaki3',
                N='cadetblue1', C='cornflowerblue', Q='cyan', S='blue', T='darkslateblue',
                R='green', H='green4', K='seagreen1',
                D='purple', E='pink',
                G='antiquewhite2', P='black', X='grey')

# Generic ggplot theme - clean with centered title by default
theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank()))

# Function to format study names
format_study <- function(x, max_width=60){
  spl <- str_split(x, '[_.]')[[1]]
  
  # Year
  year <- spl[2]
  
  # Authour
  if (all(spl[1:2] == c('roscoe', '2014'))){
    authour <- 'Roscoe & Bolon'
  } else {
    authour <- str_c(str_to_title(spl[1]), ' et al.')
  }
  
  # Gene
  gene <- spl[3]
  if (gene %in% c('nucleoprotein', 'proteing')){
    gene <- c(nucleoprotein='NP', proteing='ProteinG')[gene]
  } else if (gene == 'aph3'){
    gene <- 'APH3-II'
    spl <- spl[c(1, 2, 3, 5:length(spl))] # remove erroneously split ii
  } else {
    gene <- str_to_upper(gene)
  }
  
  # Condition
  condition <- NA
  if (!is.na(spl[4])){
    condition <- str_c(spl[4:length(spl)], collapse = ' ')
  }
  
  result <- str_c(authour, ' ', year, ' (', gene, ifelse(is.na(condition), '', str_c(', ', condition)), ')') %>%
    str_wrap(width = max_width)
  
  return(result)
}
