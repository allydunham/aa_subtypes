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

lheight <- unit(20, 'mm')
lwidth <- unit(4, 'mm')
ltitlewidth <- unit(25, 'mm')

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
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal(clip = 'off') +
  labs(x = 'UMAP1', y = 'UMAP2') +
  theme(axis.title.x = element_text(hjust = 0.525),
        axis.title.y = element_text(hjust = 0.24)) +
  annotation_raster(readPNG('figures/4_figures/position_examples/adrb2_domains.png'), interpolate = TRUE, xmin=-6, xmax=-2, ymin=4, ymax=8) +
  annotation_raster(readPNG('figures/4_figures/position_examples/cxcr4_domains.png'), interpolate = TRUE, xmin=-2.5, xmax=2, ymin=3, ymax=7.5) +
  annotation_raster(readPNG('figures/4_figures/position_examples/ccr5_domains.png'), interpolate = TRUE, xmin=2.25, xmax=6.75, ymin=4.25, ymax=8.75) +
  annotate('text', x = c(-4, 0, 4), y = 8.5, label = c('ADRB2', 'CXCR4', 'CCR5'), hjust = 0.5) + 
  annotate('point', x = c(-4, 0, 4), y = 8, shape = c(15, 18, 17), size = 2) +
  annotate('line', x = c(-5.75, 6), y = 6.5, linetype = 'dashed', colour = 'grey') + 
  annotate('line', x = c(-5.75, 6), y = 5.25, linetype = 'dashed', colour = 'grey') +
  annotate('text', x = -6, y = c(6.75, 5.875, 5),
           label = c('Extracellular', 'Transmembrane', 'Cytoplasmic'),
           colour = c('#d95f02', '#7570b3', '#1b9e77'), hjust = 1)

### Panel 2 - Mean SIFT ###
p_sift <- drop_na(dms, mean_sift) %>%
  ggplot(aes(x = umap1, y = umap2, colour = mean_sift)) +
  geom_point(colour = 'grey90', shape = 20) +
  geom_point(shape = 20) +
  scale_color_distiller(type = 'seq', palette = 'RdPu', limits = c(min(dms$mean_sift), 0),
                        breaks = 0:-4) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'log<sub>10</sub>SIFT', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth))

### Panel 3 - AA hydrophobicity ### 
p_hydrophobicity <- drop_na(dms, hydrophobicity) %>%
  ggplot(aes(x = umap1, y = umap2, colour = hydrophobicity)) +
  geom_point(colour = 'grey90', shape = 20) +
  geom_point(shape = 20) +
  scale_colour_gradientn(colours = c('#4575b4', '#e0f3f8', '#fee090', '#fc8d59', '#d73027'),
                         values = rescale01(c(-0.4, 0, 0.4, 0.8, 1.2)), limits = c(-0.4, 1.201)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'Hydrophobicity', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth))

### Panel 4 - Surface Accessibility ###
p_surface_accessibility <- drop_na(dms, side_chain_rel) %>%
  ggplot(aes(x = umap1, y = umap2, colour = all_atom_abs)) +
  geom_point(colour = 'grey90', shape = 20) +
  geom_point(shape = 20) +
  scale_colour_gradientn(colours = c('#1a2a6c', '#b21f1f', '#fdbb2d'),
                         values = c(0, 0.2, 1)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = str_wrap('Surface Accessibility', 10), barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth))

### Panel 5 - Sidechain Entropy ###
p_side_entropy <- drop_na(dms, entropy_sidechain) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(entropy_sidechain, 1.5, -1.5))) +
  geom_point(colour = 'grey90', shape = 20) +
  geom_point(shape = 20) +
  scale_colour_distiller(type = 'div', palette = 'PuOr', limits = c(-1.5, 1.5)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Sidechain Entropy (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth))

### Panel 6 - Van der Waals Clash ###
p_vdw_clash <- drop_na(dms, van_der_waals_clashes) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(van_der_waals_clashes, 5, -5))) +
  geom_point(colour = 'grey90', shape = 20) +
  geom_point(shape = 20) +
  scale_colour_distiller(type = 'div', palette = 'RdYlGn', limits = c(-5, 5)) +
  scale_y_continuous(breaks = umap2_breaks) +
  coord_equal() +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Van der Waals Clashes (kj&nbsp;mol<sup>-1</sup>)', barheight = lheight, barwidth = lwidth)) +
  theme(legend.title = element_textbox_simple(minwidth = ltitlewidth, maxwidth = ltitlewidth))

### Assemble Figure ###
size <- theme(text = element_text(size = 8))
p1 <- p_transmembrane + labs(tag = 'A') + size
p2 <- p_sift + labs(tag = 'B') + size
p3 <- p_hydrophobicity + labs(tag = 'C') + size
p4 <- p_surface_accessibility + labs(tag = 'D') + size
p5 <- p_side_entropy + labs(tag = 'E') + size
p6 <- p_vdw_clash + labs(tag = 'F') + size

figure2 <- multi_panel_figure(width = 340, height = 200, columns = 6, rows = 5,
                              panel_label_type = 'none', row_spacing = 0.1, column_spacing = 0.1) %>%
  fill_panel(p1, row = 1:5, column = 1:4) %>%
  fill_panel(p2, row = 1, column = 5:6) %>%
  fill_panel(p3, row = 2, column = 5:6) %>%
  fill_panel(p4, row = 3, column = 5:6) %>%
  fill_panel(p5, row = 4, column = 5:6) %>%
  fill_panel(p6, row = 5, column = 5:6)
ggsave('figures/4_figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')
ggsave('figures/4_figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')

