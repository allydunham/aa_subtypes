#!/usr/bin/env Rscript
# Calculate the correlation between study scores and SIFT results
source('src/config.R')

sift_dir <- 'data/sift/'
study_dirs <- dir('data/studies', full.names = TRUE)

# Import a study
import_study <- function(d){
  study <- str_split(d, '/')[[1]]
  study <- study[length(study)]
  yaml = read_yaml(str_c(d, '/', study, '.yaml'))
  
  tbl <- read_tsv(str_c(d, '/', study, '.tsv')) %>%
    mutate(study = yaml$study,
           gene = yaml$gene)
  return(tbl)
}

# Import sift results
import_sift <- function(gene){
  gene <- gene_to_filename(gene)
  fa <- as.character(readAAStringSet(str_c(sift_dir, '/', gene, '.fa'), format = 'fasta')[[1]])
  sift <- read_table(str_c(sift_dir, '/', gene, '.SIFTprediction'), skip = 5, comment = '//',
                     col_names = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                                   'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
                                   'V', 'W', 'X', 'Y', 'Z', '*', '-'),
                     col_types = cols(.default = col_double())) %>%
    pivot_longer(everything(), names_to = 'mut', values_to = 'sift') %>%
    mutate(position = rep(1:nchar(fa), each = 25),
           wt = str_split(fa, '')[[1]][position],
           log10_sift = log10(sift + 0.00005)) # SIFT goes to 4dp so 0.00005 is smaller than everything else
  return(sift)
}

dms <- lapply(study_dirs, import_study) %>%
  bind_rows()

sift <- sapply(unique(dms$gene), import_sift, simplify = FALSE) %>%
  bind_rows(.id = 'gene') 

dms <- left_join(dms, sift, by = c('gene', 'position', 'wt', 'mut'))

sift_correlations <- bind_rows(group_by(dms, study) %>% 
                                 do(tidy(cor.test(.$score, .$log10_sift, method = 'kendall'))),
                               group_by(dms, study) %>% 
                                 do(tidy(cor.test(.$score, .$log10_sift, method = 'pearson')))) %>%
  mutate(study_pretty = sapply(study, format_study, USE.NAMES = FALSE),
         p_cat = pretty_p_values(p.value, breaks = c(1e-48, 1e-12, 1e-06, 1e-3, 0.01, 0.05)))
  
p_sift_cor <- ggplot(sift_correlations, aes(x = study_pretty, y = estimate, fill = p_cat)) +
  facet_wrap(~method, ncol = 1) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0) +
  ggtitle('Correlation between Normalised Score and log10(SIFT)') +
  xlab('') +
  ylab('Correlation') +
  scale_fill_viridis_d(guide=guide_legend(title='p-value'), drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('figures/0_data_properties/sift_score_correlation.pdf', p_sift_cor, width = 20, height = 20, units = 'cm')
