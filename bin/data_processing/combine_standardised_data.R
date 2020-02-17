#!/usr/bin/env Rscript
# Combine all standardised datasets
source('src/config.R')
source('src/study_standardising.R')

# study_dirs <- dir('data/studies', full.names = TRUE) # For interactive use
study_dirs <- commandArgs(trailingOnly = TRUE)

dms <- sapply(study_dirs, import_study, fields = c('gene'), simplify = FALSE, filter=TRUE) %>%
  bind_rows() %>%
  group_by(study, gene, position, wt) %>%
  filter(sum(!mut == wt) >= 15) %>% # Only keep positions with a maximum of 4 missing scores
  ungroup()

# Import Sift results
sift <- sapply(unique(dms$gene), import_sift, simplify = FALSE) %>%
  bind_rows(.id = 'gene') 

# Import FoldX results
structure_config <- read_yaml('meta/structures.yaml')
import_fx_gene <- function(x){
  x <- gene_to_filename(x)
  import_foldx(str_c('data/foldx/', x, '/', 'average_', x, '.fxout'),
               structure_config[[x]]$sections)
}
foldx <- sapply(unique(dms$gene), import_fx_gene, simplify = FALSE) %>%
  bind_rows(.id = 'gene')

# Import backbone angles
angles <- sapply(unique(dms$gene),
                 function(gene){read_tsv(str_c('data/backbone_angles/', gene_to_filename(gene), '.tsv'), comment = '#')},
                 simplify = FALSE) %>%
  bind_rows(.id = 'gene') %>%
  select(gene, wt=aa, position, phi, psi)

# Import naccess data
import_nacc_gene <- function(x){
  x <- gene_to_filename(x)
  import_naccess(str_c('data/surface_accessibility/', x, '.rsa'),
                 structure_config[[x]]$sections)
}
surface_accessibility <- sapply(unique(dms$gene), import_nacc_gene, simplify = FALSE) %>%
  bind_rows(.id = 'gene') %>%
  select(-chain)

# Import chemical environment profiles
import_chem_env_gene <- function(x, chem_env='within_10.0'){
  x <- gene_to_filename(x)
  import_chem_env(str_c('data/chemical_environment/', x, '_', chem_env, '.tsv'), structure_config[[x]]$sections)
}
chemical_environments <- full_join(
  sapply(unique(dms$gene), import_chem_env_gene, chem_env='within_10.0', simplify = FALSE) %>%
    bind_rows(.id = 'gene') %>%
    select(-chain),
  sapply(unique(dms$gene), import_chem_env_gene, chem_env='aa_distance', simplify = FALSE) %>%
    bind_rows(.id = 'gene') %>%
    select(-chain)
)

# Import residue hydrophpbicity
hydrophobicity <- read_tsv('meta/residue_hydrophobicity.tsv', comment = '#') %>%
  rename(wt=AA, hydrophobicity = TW) %>% # use scale from Bandyopadhyay & Mehler 2008
  select(wt, hydrophobicity)

# Import Porter5 results
import_porter5_gene <- function(gene){
  gene <- gene_to_filename(gene)
  import_porter5(str_c('data/porter5/', gene, '.ss8'))
}
porter5 <- sapply(unique(dms$gene), import_porter5_gene, simplify = FALSE) %>%
  bind_rows(.id = 'gene') %>%
  select(-ss)

# Filter incomplete positions and impute
imputed_dms <- select(dms, -transformed_score, -raw_score, -class) %>%
  pivot_wider(id_cols = c(study, gene, position, wt, mut), names_from = mut, values_from = score) %>%
  pivot_longer(c(-study, -position, -wt, -gene), names_to = 'mut', values_to = 'score') %>%
  group_by(wt, mut) %>%
  mutate(score = ifelse(is.na(score), ifelse(wt == mut, 0, median(score, na.rm = TRUE)), score)) %>%
  ungroup()

# Combine data
dms <- rename(imputed_dms, imputed_score=score) %>%
  left_join(dms, by = c('study', 'position', 'wt', 'mut', 'gene')) %>%
  left_join(sift, by = c('gene', 'position', 'wt', 'mut')) %>%
  left_join(foldx, by = c('gene', 'position', 'wt', 'mut')) %>%
  left_join(angles, by = c('gene', 'position', 'wt')) %>%
  left_join(surface_accessibility, by = c('gene', 'position', 'wt')) %>%
  left_join(chemical_environments, by = c('gene', 'position', 'wt')) %>%
  left_join(porter5, by = c('gene', 'position', 'wt')) %>%
  left_join(hydrophobicity, by = c('wt')) %>%
  mutate(class = get_variant_class(wt, mut)) %>%
  arrange(study, position, wt, mut)

dms_wide <- make_dms_wide(dms)
pca <- tibble_pca(dms_wide, A:Y)
dms_wide <- bind_cols(dms_wide, as_tibble(pca$x))
tsne <- tibble_tsne(dms_wide, A:Y)
dms_wide <- tsne$tbl
dms_umap <- tibble_to_matrix(dms_wide, A:Y) %>%
  umap(metric = 'manhattan')
dms_wide <- bind_cols(dms_wide, set_colnames(dms_umap, c('umap1', 'umap2')) %>% as_tibble())

# Write combined data
write_tsv(dms, 'data/long_combined_mutational_scans.tsv')
write_tsv(dms_wide, 'data/combined_mutational_scans.tsv')
