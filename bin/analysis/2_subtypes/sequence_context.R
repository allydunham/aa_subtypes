#!/usr/bin/env Rscript
# Analyse subtype sequence context
library(Biostrings)
source('src/config.R')
library(ggseqlogo)

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

fasta <- map(dir('data/fasta', full.names = TRUE), readAAStringSet) %>% reduce(c)

extract_seq_context <- function(seq, window=10){
  s <- as.matrix(seq)[,1]
  s <- c(rep('-', window), s, rep('-', window))
  w <- map_chr((window + 1):(length(s) - window), ~str_c(s[(. - window):(. + window)], collapse = ''))
  tibble(position = 1:length(seq), wt = as.matrix(seq)[,1], seq_context = w)
}

windows <- lapply(fasta, extract_seq_context) %>%
  bind_rows(.id = 'geneid')

build_profiles <- function(tbl, ...){
  seq <- AAStringSet(tbl$seq_context)
  consensusMatrix(seq)
}

cluster_contexts <- select(dms, cluster, study, gene, position, wt) %>%
  mutate(geneid = gene_to_filename(gene)) %>%
  left_join(windows, by = c('geneid', 'position', 'wt')) %>%
  select(-geneid) %>%
  group_by(cluster)
cluster_contexts <- group_map(cluster_contexts, build_profiles) %>%
  set_names(group_keys(cluster_contexts)$cluster)

context_plots <- map(sort(unique(dms$wt)),
                     ~labeled_plot(
                       ggseqlogo(cluster_contexts[str_starts(names(cluster_contexts), .)], method = 'probability'),
                       unit = 'cm', height = 20, width = 30)
                     ) %>%
  set_names(sort(unique(dms$wt)))

save_plotlist(context_plots, 'figures/2_subtypes/final_subtypes/sequence_contexts')
