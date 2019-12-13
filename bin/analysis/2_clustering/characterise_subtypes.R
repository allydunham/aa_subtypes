#!/usr/bin/env Rscript
# Characterise generated clusters
source('src/config.R')
source('src/clustering.R')
library(argparser)

### Parse args and setup ###
parser <- arg_parser(description = 'Characterise AA subtypes', name = 'AA Subtype Characerisation')
parser <- add_argument(parser, arg = 'subtypes', help = 'TSV file assigning positions to clusters')
parser <- add_argument(parser, arg = '--dms', help = 'Path to DMS data', default = 'data/combined_mutational_scans.tsv')
parser <- add_argument(parser, arg = '--figures', help = 'Root directory to save figures', default = 'figures/2_clustering')
args <- parse_args(parser)

subtypes <- read_tsv(args$subtypes)
dms <- read_tsv(args$dms) %>%
  left_join(subtypes, ., by = c('study', 'gene', 'position', 'wt'))


