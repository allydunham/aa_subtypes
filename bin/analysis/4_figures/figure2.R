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

### Panel 1 - AA hydrophobicity ### 
p_hydrophobicity <- drop_na(dms, hydrophobicity) %>%
  ggplot(aes(x = umap1, y = umap2, colour = hydrophobicity)) +
  geom_point() +
  scale_colour_gradientn(colours = c('#4575b4', '#e0f3f8', '#fee090', '#fc8d59', '#d73027'),
                         values = rescale01(c(-0.4, 0, 0.4, 0.8, 1.2))) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = 'Hydrophobicity'))

### Panel 2 - Surface Accessibility ###
p_surface_accessibility <- drop_na(dms, side_chain_rel) %>%
  ggplot(aes(x = umap1, y = umap2, colour = all_atom_abs)) +
  geom_point() +
  scale_colour_gradientn(colours = c('#1a2a6c', '#b21f1f', '#fdbb2d'),
                         values = c(0, 0.2, 1)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colourbar(title = str_wrap('Surface Accessibility', 10)))

### Panel 4 - Sidechain Entropy ###
p_side_entropy <- drop_na(dms, entropy_sidechain) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(entropy_sidechain, 1.5, -1.5))) +
  geom_point() +
  scale_colour_distiller(type = 'div', palette = 'PuOr', limits = c(-1.5, 1.5)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Sidechain<br>Entropy<br>(kj mol<sup>-1</sup>)')) +
  theme(legend.title = element_markdown())

### Panel 3 - Sidechain Entropy ###
p_vdw_clash <- drop_na(dms, van_der_waals_clashes) %>%
  ggplot(aes(x = umap1, y = umap2, colour = clamp(van_der_waals_clashes, 5, -5))) +
  geom_point() +
  scale_colour_distiller(type = 'div', palette = 'RdYlGn', limits = c(-5, 5)) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_colorbar(title = 'Sidechain<br>Entropy<br>(kj mol<sup>-1</sup>)')) +
  theme(legend.title = element_markdown())

### Panel 5 - Transmembrane Domains ###
dms_domains <- left_join(dms, select(domains, uniprot_id, start, end, domain=name), by = 'uniprot_id') %>%
  filter(position <= end, position >= start) %>%
  select(gene, position, wt, domain) %>%
  distinct() %>%
  left_join(dms, ., by = c('gene', 'position', 'wt'))

p_transmembrane <- ggplot(mapping = aes(x=umap1, y=umap2)) + 
  geom_point(data = dms, colour = 'grey90', shape = 20) +
  geom_point(data = filter(dms_domains, gene %in% c('ADRB2', 'CCR5', 'CXCR4')), mapping = aes(shape = gene, colour = domain)) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  scale_shape_manual(values = c(ADRB2 = 15, CCR5 = 17, CXCR4 = 18)) + 
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_legend(title = 'Domain'), shape = guide_legend(title = 'Gene'))

### Panel 6 - Other Domains ###


### Assemble Figure ###
size <- theme(text = element_text(size = 10))
p1 <- p_mean_score + labs(tag = 'A') + size
p2 <- p_surface_accessibility + labs(tag = 'B') + size
p3 <- p_hydrophobicity + labs(tag = 'C') + size
p4 <- p_domains + labs(tag = 'D') + size

figure2 <- multi_panel_figure(width = 300, height = 200, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0.1) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2)
ggsave('figures/4_figures/figure2.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')
ggsave('figures/4_figures/figure2.png', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')

