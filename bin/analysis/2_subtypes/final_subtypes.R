#!/usr/bin/env Rscript
# Analyse the final set of chosen subtypes
source('src/config.R')
source('src/subtype_characterisation.R')
plots <- list()

sections <- map(read_yaml('meta/structures.yaml'), extract2, 'sections')
dms <- left_join(rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0.tsv'), cluster_ds0 = cluster),
                 rename(read_tsv('data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1.tsv'), cluster_ds1 = cluster),
                 by = c("study", "gene", "position", "wt")) %>%
  left_join(read_tsv('data/combined_mutational_scans.tsv'), by = c("study", "gene", "position", "wt")) %>%
  select(cluster_ds0, cluster_ds1, everything())

pdb_pos <- dms_pdb_positions(dms, sections)
dms <- mutate(dms, pdb_position = pdb_pos$position, pdb_chain = pdb_pos$chain)

full_characterisation <- full_cluster_characterisation(select(dms, cluster = cluster_ds0, everything()))
