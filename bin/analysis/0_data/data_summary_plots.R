#!/usr/bin/env Rscript
# Make summary plots of the project data
source('src/config.R')
source('src/study_standardising.R')

#### Studies ####
study_summary <- read_tsv('meta/study_summary.tsv') %>%
  rename_all(.funs = ~str_to_lower(.) %>% str_replace_all(' ', '_') %>% str_replace('_\\(missense\\)', '')) %>%
  mutate(study_pretty = sapply(study, format_study, USE.NAMES = FALSE))

filtered <- structure(ifelse(study_summary$filtered, 'red', 'black'), names = study_summary$study_pretty)

p_studies <- select(study_summary, study, study_pretty, completeness, coverage, mutated_positions) %>%
  pivot_longer(one_of('completeness', 'coverage'), names_to = 'metric') %>%
  mutate(study_pretty = add_markdown(study_pretty, colour = filtered)) %>%
  ggplot(aes(x = study_pretty, y = value, fill = mutated_positions)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, labeller = as_labeller(c(coverage='Proportion of Missense Variants', completeness='Proportion of Positions'))) +
  labs(x='', y='', title = 'Study Variants Summary') +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = guide_colourbar(title = '# Positions')) +
  scale_fill_viridis_c()
ggsave('figures/0_data/study_variants_summary.pdf', p_studies, units = 'cm', width = 20, height = 20)
########

#### Genes ####
gene_summary <- read_tsv('meta/gene_summary.tsv') %>%
  rename_all(.funs = ~str_to_lower(.) %>% str_replace_all(' ', '_') %>% str_replace('_\\(missense\\)', ''))

p_genes <- select(gene_summary, gene, completeness, coverage, mutated_positions) %>%
  pivot_longer(one_of('completeness', 'coverage'), names_to = 'metric') %>%
  ggplot(aes(x = gene, y = value, fill = mutated_positions)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, labeller = as_labeller(c(coverage='Proportion of Missense Variants', completeness='Proportion of Positions'))) +
  labs(x='', y='', title = 'Gene Variants Summary') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = guide_colourbar(title = '# Positions')) +
  scale_fill_viridis_c()
ggsave('figures/0_data/gene_variants_summary.pdf', p_genes, units = 'cm', width = 20, height = 20)
########

#### All position coverage ####
studies <- sapply(dir('data/studies', full.names = TRUE), import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows()

study_coverage <- group_by(studies, study, position, wt) %>%
  summarise(coverage = sum(!mut == '*' & !mut == wt)/19,
            synonymous = first(wt) %in% mut,
            nonsense = '*' %in% mut)

p_coverage <- ggplot(study_coverage, aes(x = position, y = coverage, fill = wt)) +
  facet_wrap(~study, labeller = as_labeller(sapply(unique(study_coverage$study), format_study, mark_filtered=TRUE)),
             scales = 'free_x', nrow = 5) +
  geom_col() +
  lims(y = c(0, 1)) +
  scale_fill_manual(values = AA_COLOURS) +
  theme(legend.position = 'bottom')
ggsave('figures/0_data/position_coverage.pdf', p_coverage, units = 'cm', width = 50, height = 30)

########

