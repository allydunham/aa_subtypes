#!/usr/bin/env Rscript
# Produce figure 2 (Mutational Landscape)
source('src/config.R')

domains <- read_tsv('meta/uniprot_domains.gff', comment = '#', 
                    col_names = c('uniprot_id', 'source', 'feature', 'start', 'end', 'tmp1', 'tmp2', 'tmp3', 'info', 'tmp4')) %>%
  select(-starts_with('tmp')) %>%
  filter(feature %in% c('Domain', 'Transmembrane', 'Topological domain')) %>%
  mutate(name = str_extract(info, 'Note=[^;]*;?') %>% str_remove(';$') %>% str_remove(' [0-9]*$') %>% str_remove('^Note=') %>% str_remove('%.*$'),
         name = ifelse(name == 'Helical', 'Transmembrane', name))

dms <- read_tsv('data/combined_mutational_scans.tsv') %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))

umap2_breaks <- c(-2.5, 0, 2.5)

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

p_transmembrane <- ggplot(mapping = aes(x=umap1, y=umap2)) + 
  geom_point(data = dms, colour = 'grey90', shape = 20) +
  geom_point(data = filter(dms_domains, gene %in% c('ADRB2', 'CCR5', 'CXCR4')), mapping = aes(shape = gene, colour = domain), show.legend = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  scale_shape_manual(values = c(ADRB2 = 15, CCR5 = 17, CXCR4 = 18)) + 
  scale_y_continuous(breaks = umap2_breaks, limits = c(-2.7, 7.9)) +
  scale_x_continuous(limits = c(-6.5, 6.5)) +
  coord_equal(ratio = 1.5, clip = 'off') +
  labs(x = 'UMAP1', y = 'UMAP2') +
  theme(axis.title.x = element_text(hjust = 0.525),
        axis.title.y = element_text(hjust = 0.24)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/adrb2_domains.png'), interpolate = TRUE, xmin=-6, xmax=-1.05, ymin=4, ymax=7.3) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cxcr4_domains.png'), interpolate = TRUE, xmin=-2.6, xmax=2.5, ymin=3.4, ymax=6.8) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_domains.png'), interpolate = TRUE, xmin=2, xmax=7.1, ymin=4.2, ymax=7.65) +
  annotate('richtext', x = c(-3.6, 0, 4.3), y = 7.5, label = c('ADRB2', 'CXCR4', 'CCR5'), hjust = 0.5, fill = NA, label.color = NA, size = 3.3) + 
  annotate('point', x = c(-3.6, 0, 4.3), y = 7, shape = c(15, 18, 17), size = 2) +
  annotate('line', x = c(-5.25, 6), y = 6, linetype = 'dashed', colour = 'grey') + 
  annotate('line', x = c(-5.25, 6), y = 5.1, linetype = 'dashed', colour = 'grey') +
  annotate('richtext', x = -5.5, y = c(6.25, 5.5, 4.75), hjust = 1, fill = NA, label.color = NA,
           label = c('Extracellular', 'Transmembrane', 'Cytoplasmic'), size = 3.1,
           colour = c('#d95f02', '#7570b3', '#1b9e77'))

### Panel 2 - Mean SIFT ###
p_sift <- drop_na(dms, mean_sift) %>%
  ggplot(aes(x = umap1, y = umap2, colour = mean_sift)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_color_distiller(type = 'seq', palette = 'RdPu', limits = c(min(dms$mean_sift), 0),
                        breaks = 0:-4) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'log<sub>10</sub>SIFT', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 3 - AA hydrophobicity ### 
p_hydrophobicity <- drop_na(dms, hydrophobicity) %>%
  ggplot(aes(x = umap1, y = umap2, colour = hydrophobicity)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_gradientn(colours = c('#4575b4', '#e0f3f8', '#fee090', '#fc8d59', '#d73027'),
                         values = rescale01(c(-0.4, 0, 0.4, 0.8, 1.2)), limits = c(-0.4, 1.201)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'Hydrophobicity', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 4 - Surface Accessibility ###
p_surface_accessibility <- drop_na(dms, side_chain_rel) %>%
  ggplot(aes(x = umap1, y = umap2, colour = all_atom_abs)) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_gradientn(colours = c('#1a2a6c', '#b21f1f', '#fdbb2d'),
                         values = c(0, 0.2, 1)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = str_wrap('Surface Accessibility', 10), barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 5 - Sidechain Entropy ###
p_side_entropy <- drop_na(dms, entropy_sidechain) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(entropy_sidechain, 1.5, -1.5))) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_distiller(type = 'div', palette = 'PuOr', limits = c(-1.5, 1.5)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Sidechain Entropy (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Panel 6 - Van der Waals Clash ###
p_vdw_clash <- drop_na(dms, van_der_waals_clashes) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(van_der_waals_clashes, 5, -5))) +
  geom_point(colour = 'grey90', shape = 20, size = 0.8) +
  geom_point(shape = 20, size = 0.8) +
  scale_colour_distiller(type = 'div', palette = 'RdYlGn', limits = c(-5, 5)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Van der Waals Clashes (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth, size = ltitlesize))

### Assemble Figure ###
size <- theme(text = element_text(size = 7))
p1 <- p_transmembrane + labs(tag = 'A') + size
p2 <- p_sift + guides(colour = FALSE) + labs(tag = 'B') + size
p2_legend <- get_legend(p_sift + size) %>% as_ggplot()
p3 <- p_hydrophobicity + guides(colour = FALSE) + labs(tag = 'C') + size
p3_legend <- get_legend(p_hydrophobicity + size) %>% as_ggplot()
p4 <- p_surface_accessibility + guides(colour = FALSE) + labs(tag = 'D') + size
p4_legend <- get_legend(p_surface_accessibility + size) %>% as_ggplot()
p5 <- p_side_entropy + guides(colour = FALSE) + labs(tag = 'E') + size
p5_legend <- get_legend(p_side_entropy + size) %>% as_ggplot()
p6 <- p_vdw_clash + guides(colour = FALSE) + labs(tag = 'F') + size
p6_legend <- get_legend(p_vdw_clash + size) %>% as_ggplot()

figure2 <- multi_panel_figure(width = c(125, 43, 15), height = 150, rows = 5,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1:5, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p2_legend, row = 1, column = 3) %>%
  fill_panel(p3, row = 2, column = 2) %>%
  fill_panel(p3_legend, row = 2, column = 3) %>%
  fill_panel(p4, row = 3, column = 2) %>%
  fill_panel(p4_legend, row = 3, column = 3) %>%
  fill_panel(p5, row = 4, column = 2) %>%
  fill_panel(p5_legend, row = 4, column = 3) %>%
  fill_panel(p6, row = 5, column = 2) %>%
  fill_panel(p6_legend, row = 5, column = 3)
ggsave('figures/4_figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')
ggsave('figures/4_figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')

