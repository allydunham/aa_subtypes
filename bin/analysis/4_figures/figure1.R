#!/usr/bin/env Rscript
# Produce figure 1 (Dataset summary)
source('src/config.R')
source('src/study_standardising.R')
data("BLOSUM62")

blosum62 <- as_tibble(BLOSUM62, rownames = 'wt') %>%
  pivot_longer(-wt, names_to = 'mut', values_to = 'blosum62')

raw <- sapply(dir('data/studies/', full.names = TRUE), import_study, fields = c('gene'), simplify = FALSE, filter=TRUE) %>%
  bind_rows() %>%
  group_by(study, gene, position, wt) %>%
  filter(sum(!mut == wt) >= 15) %>% # Only keep positions with a maximum of 4 missing scores
  ungroup()

dms <- read_tsv('data/combined_mutational_scans.tsv')

### Panel 1 - Gene Summary ###
gene_summary <- group_by(dms, gene) %>%
  summarise(n = n_distinct(position),
            n_struct = n_distinct(position[!is.na(total_energy)])) %>%
  mutate(prop = str_c(n_struct, ' / ', n),
         gene = as.factor(gene))

gene_labs <- levels(gene_summary$gene)
prop_labs <- structure(gene_summary$prop, names = as.character(gene_summary$gene))[gene_labs]

p_genes <- ggplot(gene_summary, aes(x = as.integer(gene))) +
  geom_col(aes(y = n), fill = 'cornflowerblue') +
  geom_errorbar(aes(ymin = n_struct, ymax = n_struct), colour='white', width = 0.5) +
  scale_x_continuous(breaks = 1:length(gene_labs), labels = gene_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(prop_labs),
                                         labels = prop_labs,
                                         name = 'With structures / Total positions')) + 
  coord_flip(clip = 'off') +
  labs(y = 'Number of Positions') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
        axis.ticks.y = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.y.right = element_text(margin = margin(l = 10)))

### Panel 2 - Normalisation Procedure ###
p_norm <- ggplot(tibble(x = 1, y = 1), aes(x=x, y=y)) +
  geom_point(colour = 'white') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank()) +
  annotate('text', x = 1, y = 1, label = 'Normalisation Procedure')

### Panel 3 - Replicate Studies ###
tem1_studys <- sapply(c('data/studies/firnberg_2014_tem1', 'data/studies/steinberg_2016_tem1'),
                      import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_rep_tem1 <- ggplot(tem1_studys, aes(x=firnberg_2014_tem1, y=steinberg_2016_tem1, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Firnberg et al. 2014', y = 'Steinberg & Ostermeier 2016', title = 'TEM1') +
  guides(colour = guide_legend(title = 'Variant Type'))

ubi_studys <- sapply(c('data/studies/roscoe_2013_ubi', 'data/studies/roscoe_2014_ubi'),
                     import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_rep_ubi <- ggplot(ubi_studys, aes(x=roscoe_2013_ubi, y=roscoe_2014_ubi, colour=class)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  coord_equal() +
  labs(x = 'Roscoe et al. 2013', y = 'Roscoe & Bolon 2014', title = 'UBI') +
  guides(colour = guide_legend(title = 'Variant Type'))

### Panel 4 - Blosum Correlation ###
blosum_cor <- group_by(raw, wt, mut) %>%
  summarise(score = mean(score)) %>%
  left_join(., blosum62, by = c('wt', 'mut'))

p_blosum <- ggplot(blosum_cor, aes(x = blosum62, y = score)) +
  geom_jitter(width = 0.2) +
  labs(x = 'BLOSUM62', y = 'Mean Normalised ER')

### Figure Assembly ###
size <- theme(text = element_text(size = 10))
p1 <- p_genes + labs(tag = 'A') + size
p2 <- p_norm + labs(tag = 'B') + size
p3_ubi <- p_rep_ubi + labs(tag = 'C') + size
p3_legend <- as_ggplot(get_legend(p3_ubi))
p3_ubi <- p3_ubi + guides(colour = FALSE)
p3_tem1 <- p_rep_tem1 + size + guides(colour = FALSE)
p4 <- p_blosum + labs(tag = 'D') + size

figure1 <- multi_panel_figure(width = 300, height = 200, columns = 8, rows = 2,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1:3) %>%
  fill_panel(p2, row = 1, column = 4:8) %>%
  fill_panel(p3_ubi, row = 2, column = 1:2) %>%
  fill_panel(p3_tem1, row = 2, column = 3:4) %>%
  fill_panel(p3_legend, row = 2, column = 5) %>%
  fill_panel(p4, row = 2, column = 6:8)
ggsave('figures/4_figures/figure1.pdf', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm')
ggsave('figures/4_figures/figure1.png', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm')
