#!/usr/bin/env Rscript
# figure S8 (Analyse SIFT score UMAP/PC landscape)
source('src/config.R')

domains <- read_tsv('meta/uniprot_domains.gff', comment = '#', 
                    col_names = c('uniprot_id', 'source', 'feature', 'start', 'end', 'tmp1', 'tmp2', 'tmp3', 'info', 'tmp4')) %>%
  select(-starts_with('tmp')) %>%
  filter(feature %in% c('Domain', 'Transmembrane', 'Topological domain')) %>%
  mutate(name = str_extract(info, 'Note=[^;]*;?') %>% str_remove(';$') %>% str_remove(' [0-9]*$') %>% str_remove('^Note=') %>% str_remove('%.*$'),
         name = ifelse(name == 'Helical', 'Transmembrane', name))

dms_long <- read_tsv('data/long_combined_mutational_scans.tsv')
dms <- read_tsv('data/combined_mutational_scans.tsv') %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))

sift_umap <- select(dms_long, study, gene, position, wt, mut, sift) %>%
  filter(!mut == '*') %>%
  pivot_wider(names_from = mut, values_from = sift) %>%
  tibble_to_matrix(A:Y) %>%
  umap(metric = 'manhattan')
sift_pca <- tibble_pca(dms, log10_sift_A:log10_sift_Y)

dms <- bind_cols(dms,
                 set_colnames(sift_umap, c('sift_umap1', 'sift_umap2')) %>% as_tibble(),
                 as_tibble(sift_pca$x) %>% rename_all(~str_c('sift_', .)))

lheight <- unit(14, 'mm')
lwidth <- unit(4, 'mm')
ltitlewidth <- unit(15, 'mm')
ltitlesize <- 6

### Panel 1 - Transmembrane Domains ###
dms_domains <- left_join(dms, select(domains, uniprot_id, start, end, domain=name), by = 'uniprot_id') %>%
  filter(position <= end, position >= start) %>%
  select(gene, position, wt, domain) %>%
  distinct() %>%
  left_join(dms, ., by = c('gene', 'position', 'wt'))

p_transmembrane <- ggplot(filter(dms_domains, gene %in% c('ADRB2', 'CCR5', 'CXCR4')), aes(x=sift_umap1, y=sift_umap2, shape = gene, colour = domain)) + 
  facet_wrap(~domain, nrow = 1) +
  geom_point(data = dms, colour = 'grey90', shape = 20) +
  geom_point() +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  scale_shape_manual(values = c(ADRB2 = 15, CCR5 = 17, CXCR4 = 18)) + 
  labs(x = 'UMAP1', y = 'UMAP2') +
  guides(shape = guide_legend(title = ''), colour = FALSE)

### Panel 2 - Mean Score ###
p_sift <- drop_na(dms, mean_sift) %>%
  ggplot(aes(x = sift_umap1, y = sift_umap2, colour = mean_sift)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_color_distiller(type = 'seq', palette = 'RdPu', limits = c(min(dms$mean_sift), 0),
                        breaks = 0:-4) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'log<sub>10</sub>SIFT', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 3 - AA hydrophobicity ### 
p_hydrophobicity <- drop_na(dms, hydrophobicity) %>%
  ggplot(aes(x = sift_umap1, y = sift_umap2, colour = hydrophobicity)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_gradientn(colours = c('#4575b4', '#e0f3f8', '#fee090', '#fc8d59', '#d73027'),
                         values = rescale01(c(-0.4, 0, 0.4, 0.8, 1.2)), limits = c(-0.4, 1.201)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'Hydrophobicity', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 4 - Surface Accessibility ###
p_surface_accessibility <- drop_na(dms, all_atom_abs) %>%
  ggplot(aes(x = sift_umap1, y = sift_umap2, colour = all_atom_abs)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_gradientn(colours = c('#1a2a6c', '#b21f1f', '#fdbb2d'),
                         values = c(0, 0.2, 1)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = str_wrap('Surface Accessibility', 10), barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 5 - Sidechain Entropy ###
p_side_entropy <- drop_na(dms, entropy_sidechain) %>%
  ggplot(aes(x = sift_umap1, y = sift_umap2, colour = clamp(entropy_sidechain, 1.5, -1.5))) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_distiller(type = 'div', palette = 'PuOr', limits = c(-1.5, 1.5)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Sidechain Entropy (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 6 - Van der Waals Clash ###
p_vdw_clash <- drop_na(dms, van_der_waals_clashes) %>%
  ggplot(aes(x = sift_umap1, y = sift_umap2, colour = clamp(van_der_waals_clashes, 5, -5))) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_distiller(type = 'div', palette = 'RdYlGn', limits = c(-5, 5)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Van der Waals Clashes (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 7 - Score correlation ###
p_pc_score <- ggplot(dms, aes(x = sift_PC1, y = mean_sift, z = mean_score)) +
  stat_summary_hex(bins = 40) +
  scale_fill_distiller(type = 'div', palette = 'Spectral', direction = 1) +
  labs(x = 'log<sub>10</sub>SIFT PC1', y = 'Mean log<sub>10</sub>SIFT') +
  guides(fill = guide_colourbar(title = 'Mean ER', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize),
        axis.title = element_markdown())

### Assemble figure ###
size <- theme(text = element_text(size = 7))
p1 <- p_transmembrane + labs(tag = 'A') + size
p2 <- p_sift + labs(tag = 'B') + size
p3 <- p_hydrophobicity + labs(tag = 'C') + size
p4 <- p_surface_accessibility + labs(tag = 'D') + size
p5 <- p_side_entropy + labs(tag = 'E') + size
p6 <- p_vdw_clash + labs(tag = 'F') + size
p7 <- p_pc_score + labs(tag = 'G') + size

figure <- multi_panel_figure(width = 183, height = c(55, 50, 50), columns = 3,
                             column_spacing = 0, row_spacing = 0, unit = 'mm',
                             panel_label_type = 'none') %>%
  fill_panel(p1, row = 1, column = 1:3) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p4, row = 2, column = 3) %>%
  fill_panel(p5, row = 3, column = 1) %>%
  fill_panel(p6, row = 3, column = 2) %>%
  fill_panel(p7, row = 3, column = 3)

ggsave('figures/4_figures/figureS8.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS8.png', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
