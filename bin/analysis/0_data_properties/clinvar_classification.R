#!/usr/bin/env Rscript
# Test ability of DMS data to classify clinvar pathogenic variants, compared to SIFT and FoldX scores
library(VariantAnnotation)
source('src/config.R')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg38)

# Import DMS data and identify human genes
meta <- read_tsv('meta/gene_summary.tsv')
dms <- read_tsv('data/long_combined_mutational_scans.tsv')

dms_genes <- dplyr::filter(meta, Species == 'H. sapiens') %>% 
  pull(`Uniprot ID`)

#### Import ClinVar Data ####
# Find locations of genes of interest
dms_genes_entrez_id <- select(Homo.sapiens, dms_genes, "ENTREZID", "UNIPROT") %>%
  distinct(UNIPROT, .keep_all = TRUE)

dms_gene_ranges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
dms_gene_ranges <- dms_gene_ranges[dms_gene_ranges$gene_id %in% dms_genes_entrez_id$ENTREZID]
seqlevels(dms_gene_ranges) <- structure(c(1:22, 'X', 'Y'), names=str_c('chr', c(1:22, 'X', 'Y')))

# Import clinvar data from appropriate regions
tbi <- TabixFile('data/clinvar/clinvar_20191202.vcf.gz')
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
########

#### Combine DMS and clinvar data ####
variant_data <- dplyr::select(dms, gene, position, wt, mut, score, raw_score, class, log10_sift, total_energy) %>%
  drop_na(score) %>%
  left_join(dplyr::select(meta, gene=Gene, uniprot_id=`Uniprot ID`), by = 'gene') %>%
  left_join(clinvar_variants, by = c('uniprot_id', 'position', 'wt', 'mut')) %>%
  drop_na(clnsig) %>%
  filter(!clnsig %in% c('Uncertain_significance', 'not_provided', 'drug_response', 'Conflicting_interpretations_of_pathogenicity')) %>%
  mutate(clnsig_simple = if_else(clnsig %in% c('Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic'), 'pathogenic', 'benign'),
         clnsig_train = clnsig_simple == 'pathogenic')
########

#### ROC Curves / AUC ####
calc_true_false_over_range <- function(tbl, label, preds, n=100, invert=FALSE){
  label <- enquo(label)
  preds <- enquo(preds)
  
  if (all(is.na(pull(tbl, !!preds)))){
    return(tibble(TP=NA, TN=NA, FP=NA, FN=NA))
  }
  
  r <- range(pull(tbl, !!preds), na.rm = TRUE)
  step <- (r[2] - r[1])/n
  if (invert){
    threshs <- seq(r[2], r[1], -step)
  } else {
    threshs <- seq(r[1], r[2], step)
  }
  
  
  bind_rows(sapply(threshs, function(x){calc_true_false(tbl, !!label, !!preds, x, invert)}, simplify = FALSE)) %>%
    mutate(threshold = threshs) %>%
    return()
}

calc_true_false <- function(tbl, label, preds, threshold, invert=FALSE){
  label <- enquo(label)
  preds <- enquo(preds)
  
  if (invert) {
    tbl <- dplyr::select(tbl, !!label, !!preds) %>%
      drop_na() %>%
      mutate(pred = !!preds < threshold)
  } else {
    tbl <- dplyr::select(tbl, !!label, !!preds) %>%
      drop_na() %>%
      mutate(pred = !!preds > threshold)
  }
  
  mutate(tbl,
         TP = pred & !!label,
         TN = !pred & !(!!label),
         FP = pred & !(!!label),
         FN = !pred & !!label) %>%
    summarise_at(.vars = vars(TP, TN, FP, FN), .funs = sum) %>%
    return()
}

trapezium_rule <- function(x, y){
  x <- x[order(x)]
  y <- y[order(y)]
  
  sum(diff(x) * (y[-length(y)] + y[-1])/2)
}

calc_rocs <- function(tbl){
  if (sum(!is.na(tbl$total_energy)) < 5) {
    return(calc_rocs_no_fx(tbl))
  }
  
  # Create Combined Model using all three scores
  model_foldx <- glm(clnsig_train ~ score + log10_sift + total_energy, data = tbl, family = binomial())
  model_sift <- glm(clnsig_train ~ score + log10_sift, data = tbl, family = binomial())
  
  tbl <- mutate(tbl, model_sift = model_sift$fitted.values, model_foldx = NA)
  tbl[as.integer(names(model_foldx$fitted.values)), 'model_foldx'] <- model_foldx$fitted.values
  
  # Cacl TP/FP
  bind_rows(ER=calc_true_false_over_range(tbl, clnsig_train, score, invert = TRUE),
            SIFT=calc_true_false_over_range(tbl, clnsig_train, log10_sift, invert = TRUE),
            FoldX=calc_true_false_over_range(tbl, clnsig_train, total_energy, invert = FALSE),
            `Model (ER, SIFT)`=calc_true_false_over_range(tbl, clnsig_train, model_sift, invert = FALSE),
            `Model (ER, SIFT, FoldX)`=calc_true_false_over_range(tbl, clnsig_train, model_foldx, invert = FALSE),
            .id = 'score') %>%
    mutate(TPR = TP / (TP + FN),
           FPR = FP / (FP + TN),
           Precision = TP / (TP + FP),
           Recall = TPR)
}

calc_rocs_no_fx <- function(tbl){
  model_sift <- glm(clnsig_train ~ score + log10_sift, data = tbl, family = binomial())
  
  tbl <- mutate(tbl, model_sift = model_sift$fitted.values, model_foldx = NA)
  
  # Cacl TP/FP
  bind_rows(ER=calc_true_false_over_range(tbl, clnsig_train, score, invert = TRUE),
            SIFT=calc_true_false_over_range(tbl, clnsig_train, log10_sift, invert = TRUE),
            `Model (ER, SIFT)`=calc_true_false_over_range(tbl, clnsig_train, model_sift, invert = FALSE),
            .id = 'score') %>%
    mutate(TPR = TP / (TP + FN),
           FPR = FP / (FP + TN),
           Precision = TP / (TP + FP),
           Recall = TPR)
}

variant_performance <- calc_rocs(variant_data)

# Calc AUC
score_auc <- group_by(variant_performance, score) %>%
  summarise(auc = trapezium_rule(c(0, FPR, 1), c(0, TPR, 1))) %>%
  mutate(score = factor(score, levels = score[order(auc, decreasing = TRUE)]))

score_cols <- RColorBrewer::brewer.pal(n = length(score_auc$score), name = 'Set2')
names(score_cols) <- levels(score_auc$score)

p_roc <- ggplot(variant_performance, aes(x=FPR, y=TPR, colour=score)) + 
  geom_path() +
  geom_abline(slope = 1, linetype = 'dotted') +
  coord_equal() +
  guides(colour=guide_legend(title = '')) +
  scale_color_manual(values = score_cols) +
  theme(text = element_text(size = 9))

p_auc <- ggplot(score_auc, aes(x = score, y = auc, fill = score)) +
  geom_col(width = 0.75) +
  guides(fill = FALSE) +
  labs(x = '', y = 'AUC') + 
  lims(y = c(0, 1)) +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank()) +
  scale_fill_manual(values = score_cols)

ggsave('figures/0_data_properties/clinvar_auc.pdf', p_auc, width = 4, height = 8, units = 'cm')
ggsave('figures/0_data_properties/clinvar_roc.pdf', p_roc, width = 12, height = 8, units = 'cm')

########

#### ROC / AUC separately between genes ####
per_gene_performance <- group_by(variant_data, gene) %>%
  mutate(n = n()) %>%
  filter(n > 40) %>%
  group_modify(~calc_rocs(.x))

per_gene_auc <- group_by(per_gene_performance, gene, score) %>%
  summarise(auc = trapezium_rule(c(0, FPR, 1), c(0, TPR, 1))) %>%
  mutate(score = factor(score, levels = levels(score_auc$score)))

p_roc_per_gene <- ggplot(per_gene_performance, aes(x=FPR, y=TPR, colour=score)) + 
  geom_path() +
  geom_abline(slope = 1, linetype = 'dotted') +
  facet_wrap(~gene) +
  coord_equal() +
  guides(colour=guide_legend(title = '')) +
  scale_color_manual(values = score_cols) +
  theme(text = element_text(size = 9))

p_auc_per_gene <- ggplot(per_gene_auc, aes(x = score, y = auc, fill = score)) +
  geom_col(width = 0.75) +
  facet_wrap(~gene) +
  guides(fill = FALSE) +
  labs(x = '', y = 'AUC') + 
  lims(y = c(0, 1)) +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank()) +
  scale_fill_manual(values = score_cols)

ggsave('figures/0_data_properties/clinvar_auc_per_gene.pdf', p_auc_per_gene, width = 8, height = 16, units = 'cm')
ggsave('figures/0_data_properties/clinvar_roc_per_gene.pdf', p_roc_per_gene, width = 24, height = 16, units = 'cm')

#######
