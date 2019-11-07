#!/usr/bin/env Rscript
# Make summary plots of the project data

source('src/config.R')

#### Studies ####
studies <- read_tsv('meta/study_summary.tsv') %>%
  rename_all(.funs = ~str_to_lower(.) %>% str_replace_all(' ', '_') %>% str_replace('_\\(missense\\)', '')) %>%
  mutate(study_pretty = sapply(study, format_study, USE.NAMES = FALSE))

filtered <- structure(ifelse(studies$filtered, 'red', 'black'), names = studies$study_pretty)

p_studies <- select(studies, study, study_pretty, completeness, coverage, mutated_positions) %>%
  pivot_longer(one_of('completeness', 'coverage'), names_to = 'metric') %>%
  ggplot(aes(x = study_pretty, y = value, fill = mutated_positions)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, labeller = as_labeller(c(coverage='Proportion of Missense Variants', completeness='Proportion of Positions'))) +
  labs(x='', y='', title = 'Study Variants Summary') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    colour = filtered[unique(studies$study_pretty)])) +
  guides(fill = guide_colourbar(title = '# Positions')) +
  scale_fill_viridis_c()
ggsave('figures/0_data_properties/study_variants_summary.pdf', p_studies, units = 'cm', width = 20, height = 20)
########

#### Genes ####
genes <- read_tsv('meta/gene_summary.tsv') %>%
  rename_all(.funs = ~str_to_lower(.) %>% str_replace_all(' ', '_') %>% str_replace('_\\(missense\\)', ''))

p_genes <- select(genes, gene, completeness, coverage, mutated_positions) %>%
  pivot_longer(one_of('completeness', 'coverage'), names_to = 'metric') %>%
  ggplot(aes(x = gene, y = value, fill = mutated_positions)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, labeller = as_labeller(c(coverage='Proportion of Missense Variants', completeness='Proportion of Positions'))) +
  labs(x='', y='', title = 'Gene Variants Summary') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = guide_colourbar(title = '# Positions')) +
  scale_fill_viridis_c()
ggsave('figures/0_data_properties/gene_variants_summary.pdf', p_genes, units = 'cm', width = 20, height = 20)
########