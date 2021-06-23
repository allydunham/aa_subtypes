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
dms_long <- read_tsv('data/long_combined_mutational_scans.tsv')
  
### Panel 1 - Gene Summary ###
gene_summary <- group_by(dms, gene) %>%
  summarise(n = n_distinct(position),
            n_struct = n_distinct(position[!is.na(total_energy)])) %>%
  mutate(percent_struct = n_struct / n * 100,
         x = str_c(n, " (", signif(percent_struct, 3), "%)"),
         img = str_c("<img src='figures/4_figures/proteins/", gene_to_filename(gene), ".png", "' width='53' />"))

p_genes <- ggplot(gene_summary, aes(x = x, y = n, label = img)) +
  geom_richtext(fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_wrap(~gene, nrow = 5, scales = 'free') +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

### Panel 2 - Normalisation Procedure ###
p_norm_raw <- filter(raw, study %in% c('steinberg_2016_tem1', 'heredia_2018_ccr5', 'matreyek_2018_tpmt')) %>%
  ggplot(aes(x = raw_score, colour = study)) +
  geom_density(size = 0.9) +
  labs(title = 'Raw') +
  guides(colour = FALSE) +
  scale_x_continuous(breaks = c(0),
                     labels = c("<span style='font-size:6pt'><span style='line-height:4'>0</span><br>Studies on<br>different scales</span>")) +
  theme(text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_markdown(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(5, 0, 5, 0), units = 'mm'))

p_norm_trans <- filter(raw, study %in% c('steinberg_2016_tem1', 'heredia_2018_ccr5', 'matreyek_2018_tpmt')) %>%
  ggplot(aes(x = transformed_score, colour = study)) +
  geom_density(size = 0.9) +
  labs(title = 'Transformed') + 
  guides(colour = FALSE) +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(breaks = c(0, -6.5), labels = c("<span style='font-size:6pt'><span style='line-height:4'>0</span><br>Neutral</span>", "<img src='figures/4_figures/parts/arrow.png' width='40' /><span style='font-size:6pt'><br>Deleterious</span>")) +
  scale_y_continuous(expand = expansion(0,0)) +
  theme(text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_markdown(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(5, 0, 5, 0), units = 'mm'))

p_norm_final <- filter(raw, study %in% c('steinberg_2016_tem1', 'heredia_2018_ccr5', 'matreyek_2018_tpmt')) %>%
  ggplot(aes(x = score, colour = study)) +
  geom_density(size = 0.9) +
  labs(title = 'Normalised') + 
  guides(colour = FALSE) +
  scale_x_continuous(limits = c(-2, 1), breaks = c(0, -1), labels = c('0', '-1\nMean of 10%\nmost deleterious')) +
  scale_y_continuous(expand = expansion(0,0)) +
  theme(text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(5, 0, 5, 0), units = 'mm'))

p_norm <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1.1), y = c(0, 2)) + 
  coord_cartesian(clip = 'off') +
  labs(title = 'Combining Deep Mutational\nScanning Studies') + 
  annotate('polygon', x = c(0.25, 0.5, 0.35, 0.5, 0, 0.15, 0), y = c(0.25, 0.5, 0.5, 2, 2, 0.5, 0.5), fill = 'lightskyblue') +
  annotate('richtext', label = "<span style='font-size:8pt'>Process double<br>mutants and<br>multiple mutations</span>", x = 0.25, y = 1.85, fill = NA, label.colour = NA) +
  annotate('richtext', label = "<span style='font-size:8pt'>Transform scores<br>to standard scale</span>", x = 0.25, y = 1.45, fill = NA, label.colour = NA) +
  annotate('richtext', label = "<span style='font-size:8pt'>Normalise against<br>mean of lowest<br>10% of scores</span>", x = 0.25, y = 1, fill = NA, label.colour = NA) +
  annotate('richtext', label = "<span style='font-size:8pt'>Filter positions<br>with scores for<br>&lt;15 substitutions</span>", x = 0.25, y = 0.5, fill = NA, label.colour = NA) +
  annotate('richtext', label = "<span style='font-size:10pt'>**Comparable<br>ER Scores**</span>", x = 0.25, y = 0.15, fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_norm_raw), xmin = 0.55, xmax = 1.1, ymin = 1.5, ymax = 2.1) +
  annotation_custom(ggplotGrob(p_norm_trans), xmin = 0.55, xmax = 1.1, ymin = 0.8, ymax = 1.5) + 
  annotation_custom(ggplotGrob(p_norm_final), xmin = 0.55, xmax = 1.1, ymin = 0, ymax = 0.7) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

### Panel 3 - Replicate Studies ###
ubi_studys <- sapply(c('data/studies/roscoe_2013_ubi', 'data/studies/roscoe_2014_ubi'),
                     import_study, fields = 'gene', simplify = FALSE) %>%
  bind_rows() %>%
  select(-transformed_score, -raw_score, -gene) %>%
  pivot_wider(names_from = study, values_from = score)

p_rep_ubi <- ggplot(ubi_studys, aes(x=roscoe_2013_ubi, y=roscoe_2014_ubi, colour=class)) +
  geom_point(shape=20) +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'black') +
  geom_abline(slope = 1, colour = 'black', linetype = 'dotted') +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  scale_x_continuous(expand = expansion(0.01)) +
  scale_y_continuous(expand = expansion(0.01)) +
  guides(colour = FALSE) +
  labs(x = 'Roscoe et al. 2013', y = 'Roscoe & Bolon 2014') +
  theme(panel.grid.major.y = element_blank())

### Panel 4 - Blosum Correlation ###
blosum_cor <- group_by(raw, wt, mut) %>%
  summarise(score = mean(score)) %>%
  left_join(., blosum62, by = c('wt', 'mut')) %>%
  mutate(mut_class = get_variant_class(wt, mut))

p_blosum <- ggplot(blosum_cor, aes(x = blosum62, y = score, colour = mut_class)) +
  geom_jitter(width = 0.2, shape = 20) +
  geom_smooth(data = filter(blosum_cor, mut_class == 'Missense'), mapping = aes(x = blosum62, y = score), colour = 'black', method = 'lm', formula = y ~ x, inherit.aes = FALSE) +
  scale_colour_manual(values = MUT_CLASS_COLOURS) +
  guides(colour = guide_legend(title = '', override.aes = list(size = 2), direction = 'horizontal')) +
  labs(x = 'BLOSUM62', y = 'Mean Normalised ER') +
  theme(legend.position = 'bottom')

p_variant_legend <- get_legend(p_blosum) %>% as_ggplot()
p_blosum <- p_blosum + guides(colour = FALSE)

### Panel 5 Sift correlation ###
format_gene <- function(gene, study){
  year <- str_split_fixed(study, fixed('_'), n = 3)[,2]
  ifelse(gene %in% c('UBI', 'HSP90'), str_c(gene, ' (', year, ')'), gene)
}

sift_correlations <- select(dms_long, study, gene, position, wt, mut, score, log10_sift) %>%
  drop_na(score, log10_sift) %>%
  group_by(study, gene) %>% 
  group_modify(~mutate(tidy(cor.test(.$score, .$log10_sift, method = 'pearson')), n = length(.$score))) %>%
  ungroup() %>%
  mutate(study_pretty = sapply(study, format_study, USE.NAMES = FALSE),
         p_cat = pretty_p_values(p.value, breaks = c(1e-48, 1e-12, 1e-06, 1e-3, 0.01, 0.05), markdown_exp = TRUE, prefix_p = TRUE),
         gene_pretty = as.factor(format_gene(gene, study))) %>%
  arrange(gene_pretty)

p_sift <- ggplot(sift_correlations, aes(x = as.integer(gene_pretty), y = estimate, fill = p_cat, label = str_c("n = ", n))) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.5, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, colour = 'grey') +
  scale_x_continuous(breaks = 1:length(levels(sift_correlations$gene_pretty)), labels = as.character(sift_correlations$gene_pretty),
                     sec.axis = dup_axis(name = "", labels = sift_correlations$n), expand = expansion(0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x='', y=expression("Pearson's"~rho)) +
  scale_fill_viridis_d(drop=FALSE) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 6),
        legend.key.size = unit(2, 'mm'),
        legend.margin = margin(-10,0,0,0, unit = "pt"),
        legend.box.spacing = unit(1, "mm"),
        legend.position = 'bottom')

### Figure Assembly ###
size <- theme(text = element_text(size = 8))
p1 <- p_genes + labs(tag = 'A') + size
p2 <- p_norm + labs(tag = 'B') + size
p3 <- p_rep_ubi + labs(tag = 'C') + size
p4 <- p_blosum + labs(tag = 'D') + size
p34_legend <- p_variant_legend + size
p5 <- p_sift + labs(tag = 'E') + size

figure1 <- multi_panel_figure(width = 183, height = 220, columns = 20, rows = 20,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1:14, column = 1:14) %>%
  fill_panel(p2, row = 1:14, column = 15:20) %>%
  fill_panel(p3, row = 15:19, column = 1:5) %>%
  fill_panel(p4, row = 15:19, column = 6:10) %>%
  fill_panel(p34_legend, row = 20, column = 1:10) %>%
  fill_panel(p5, row = 15:20, column = 11:20)
ggsave('figures/4_figures/figure1.pdf', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm')
ggsave('figures/4_figures/figure1.png', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm')
ggsave('figures/4_figures/figure1.tiff', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm')
ggsave('figures/4_figures/figure1.eps', figure1, width = figure_width(figure1), height = figure_height(figure1), units = 'mm', device=cairo_ps, fallback_resolution = 600)
