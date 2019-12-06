#!/usr/bin/env Rscript
# Test ability of DMS data to classify clinvar pathogenic variants, compared to SIFT and FoldX scores

library(VariantAnnotation)
library(ensembldb)
library(AnnotationHub)
source('src/config.R')

# Load GRCh38 Ensembl Genome Annotation
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb", 98))
ensembl_db <- ahDb[[1]]

# Import DMS data and clinvar VCF
meta <- read_tsv('meta/gene_summary.tsv')
dms <- read_tsv('data/long_combined_mutational_scans.tsv')
clinvar <- readVcf('data/clinvar/clinvar_20191202.vcf', "grch38")

# Identify human genes we have DMS data for
dms_genes <- filter(meta, Species == 'H. sapiens') %>% 
  pull(Gene) %>% 
  str_to_upper() %>%
  str_replace('RAS', 'HRAS') %>% # Switch to gene names for those using protein names/common names
  str_replace('TDP43', 'TARDBP')
