#!/usr/bin/env Rscript
# Summarise standardised dataset

source('src/config.R')

dms <- read_tsv('data/combined_mutational_scans.tsv')

p_score_dist <- ggplot(dms, aes(x=imputed_score, fill=ifelse(!is.na(score), 'Experiment', 'Imputed'))) +
  facet_wrap(~study, ncol = 4, labeller = as_labeller(sapply(unique(dms$study), format_study, max_width = 28)), scales = 'free') +
  geom_histogram(bins=30) +
  labs(x = 'Normalised ER', y='Count') +
  scale_fill_manual(values = c(Experiment='cornflowerblue', Imputed='firebrick2')) +
  guides(fill=guide_legend(title = ''))
ggsave('figures/0_data_properties/standardised_distributions.pdf', p_score_dist, units = 'cm', height = 35, width = 30)

summary_tbl <- group_by(dms, study, position, wt) %>%
  summarise(fx = sum(!is.na(total_energy)) == 19,
            sift = all(!is.na(sift)),
            phi = all(!is.na(phi)),
            psi = all(!is.na(psi)),
            sa = all(!is.na(all_atom_abs))) %>%
  ungroup() %>%
  summarise(Total = n(),
            `FoldX Results` = sum(fx),
            `SIFT Results` = sum(sift),
            Phi = sum(phi),
            Psi = sum(psi),
            `Surface Accessibility` = sum(sa)) %>%
  pivot_longer(everything(), names_to = 'metric', values_to = 'Count') %>%
  mutate(metric = factor(metric, levels = metric[order(Count, c(1, rep(0, length(metric) - 1)))])) # Second order arg ensures total is always first

p_position_summary <- ggplot(summary_tbl, aes(x = metric, y = Count, fill = metric)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = Count), nudge_y = 500) +
  coord_flip() +
  guides(fill=FALSE) +
  labs(x='', title = 'Summary of data collected after filtering') +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('figures/0_data_properties/position_data_summary.pdf', p_position_summary, units = 'cm', height = 8, width = 15)
