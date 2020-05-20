#!/usr/bin/env Rscript
# Produce figure SX (Disease variant subtypes)
library(VariantAnnotation)
source('src/config.R')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg38)

### Import DMS data ###
dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position) %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))

dms_genes <- unique(dms$uniprot_id)

### Import ClinVar Data ###
# Find locations of genes of interest
dms_genes_entrez_id <- select(Homo.sapiens, dms_genes, "ENTREZID", "UNIPROT") %>%
  distinct(UNIPROT, .keep_all = TRUE)

dms_gene_ranges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
dms_gene_ranges <- dms_gene_ranges[dms_gene_ranges$gene_id %in% dms_genes_entrez_id$ENTREZID]
seqlevels(dms_gene_ranges) <- structure(c(1:22, 'X', 'Y'), names=str_c('chr', c(1:22, 'X', 'Y')))

# Import clinvar data from appropriate regions
tbi <- TabixFile('data/clinvar/clinvar_20200106.vcf.gz')
clinvar <- readVcf(tbi, "hg38", param = dms_gene_ranges)
seqlevels(rowRanges(clinvar)) <- structure(str_c('chr', seqlevels(rowRanges(clinvar))), names=seqlevels(rowRanges(clinvar)))

# Identify coding variants
clinvar_coding <- predictCoding(clinvar, TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg38)
clinvar_coding <- clinvar_coding[(sapply(clinvar_coding$PROTEINLOC, length) == 1) & 
                                   (clinvar_coding$GENEID %in% dms_genes_entrez_id$ENTREZID) &
                                   (clinvar_coding$CONSEQUENCE %in% c('nonsynonymous', 'nonsense'))]

# Match Coding variants up to genes and clinvar significances
clinvar_variants <- tibble(query = clinvar_coding$QUERYID, entrez_id = clinvar_coding$GENEID, position = drop(clinvar_coding$PROTEINLOC),
                           wt = as.character(clinvar_coding$REFAA), mut = as.character(clinvar_coding$VARAA)) %>%
  distinct() %>%
  left_join(rename(dms_genes_entrez_id, entrez_id=ENTREZID, uniprot_id = UNIPROT), by = 'entrez_id') %>%
  mutate(clnsig = sapply(info(clinvar)$CLNSIG, extract, 1)[query])

### Combine DMS and clinvar data ###
genes <- unique(dms$uniprot_id)
genes <- genes[genes %in% unique(clinvar_variants$uniprot_id)]
dms <- dplyr::select(clinvar_variants, uniprot_id, position, wt, cv_mut=mut, clnsig) %>%
  left_join(dms, ., by=c('uniprot_id', 'position', 'wt')) %>%
  mutate(clnsig_simple = if_else(clnsig %in% c('Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic'), 'pathogenic', 'benign'))

### Analyse subtypes' clinvar variants ###
subtype_path <- filter(dms, !str_detect(cluster, CLUSTER_OUTLIER_RE)) %>%
  group_by(wt, cluster) %>%
  summarise(n = n(),
            path = sum(clnsig_simple == 'pathogenic')) %>%
  mutate(n_aa = sum(n),
         path_aa = sum(path)) %>%
  ungroup() %>%
  mutate(freq = path / n,
         freq_aa = path_aa / n_aa,
         ratio = log2(pmax(freq, 0.01)/freq_aa),
         p = pmap(.l = list(path, n, freq_aa), .f = ~binom.test(x = ..1, n = ..2, p = ..3)$p.value) %>% unlist(),
         p_adj = p.adjust(p, method = 'BH'))

ggplot(subtype_path, aes(x = ratio, y = p_adj, colour = p_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values = c(`TRUE`='red', `FALSE`='black'))
