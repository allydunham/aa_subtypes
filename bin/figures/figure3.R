#!/usr/bin/env Rscript
# Produce figure 3 (Subtype Overview)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

# Subtype Categories
most_selective_subtypes <- group_by(full_characterisation$profiles, cluster) %>%
  summarise(mean_er = mean(er)) %>%
  mutate(aa = str_sub(cluster, end = 1)) %>%
  group_by(aa) %>%
  filter(!mean_er > min(mean_er))

### Panel 1 - Schematic ###
## Subparts
p_initial_profiles <- map(c('A', 'C', 'D', 'W', 'Y'), ~filter(dms, wt == .) %>%
                            select(study, cluster, position, wt, A:Y) %>%
                            extract(1:20, 1:24) %>%
                            pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
                            mutate(er = clamp(er, 1, -1)) %>%
                            ggplot(aes(x = str_c(study, position), y = mut, fill = er)) +
                            geom_raster() +
                            coord_fixed() +
                            scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                                                 limits = c(-1, 1)) +
                            guides(fill = FALSE) +
                            theme(axis.text = element_blank(),
                                  axis.title = element_blank(),
                                  axis.ticks = element_blank(),
                                  panel.grid.major.y = element_blank(),
                                  plot.margin = unit(c(0, 0, 0, 0), 'mm'))
  ) %>% set_names(c('A', 'C', 'D', 'W', 'Y'))

p_a_profiles <- map(c('AP', 'Rest'), ~filter(dms, wt == 'A') %>%
                      mutate(cat = ifelse(cluster == 'AP', 'AP', 'Rest')) %>%
                      filter(cat == .x) %>%
                      select(study, cluster, position, wt, A:Y) %>%
                      extract(1:20, 1:24) %>%
                      pivot_longer(A:Y, names_to = 'mut', values_to = 'er') %>%
                      mutate(er = clamp(er, 1, -1)) %>%
                      drop_na(er) %>%
                      ggplot(aes(x = str_c(study, position), y = mut, fill = er)) +
                      geom_raster() +
                      coord_fixed() +
                      scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                                           limits = c(-1, 1)) +
                      guides(fill = FALSE) +
                      theme(axis.text = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks = element_blank(),
                            panel.grid.major.y = element_blank(),
                            plot.margin = unit(c(0, 0, 0, 0), 'mm'))
) %>% set_names(c('AP', 'Rest'))

schematic_profiles <- pivot_wider(full_characterisation$profiles, cluster, names_from = 'mut', values_from = 'er') %>% filter(str_starts(cluster, 'A'))
schematic_distance <- tibble_to_matrix(schematic_profiles, A:Y, row_names = 'cluster') %>% cosine_distance_matrix() %>% as.dist()
schematic_hc <- hclust(schematic_distance)
schematic_dend_data <- dendro_data(schematic_hc)
schematic_branches <- schematic_dend_data$segments
schematic_leaves <- schematic_dend_data$labels
p_schematic_dend <- ggplot() +
  geom_segment(data = schematic_branches, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_point(data = schematic_leaves, aes(x=x, y=y, colour=label), shape=20, size=1.5, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2', na.value = 'grey', direction = -1)

p_schematic_dend_cuts <- p_schematic_dend +
  annotate('line', x = c(1.6, 2.25), y = c(0.27, 0.28), size=1.2, colour='red') +
  annotate('line', x = c(3.45, 4.05), y = c(0.17, 0.18), size=1.2, colour='red') +
  annotate('line', x = c(6.2, 6.8), y = c(0.47, 0.48), size=1.2, colour='red')

## Full Figure
p_schematic <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1.1), y = c(0, 1.1)) + 
  coord_fixed(clip = 'off') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank()) +
  annotate('richtext', x = 0.5, y = 1.1, label = "<span style='font-size:9pt'>Clustering ER Profiles</span>", fill = NA, label.colour = NA) +
  
  # Initial profile heatmaps
  annotation_custom(ggplotGrob(p_initial_profiles$A), xmin = 0.14, xmax = 0.29, ymin = 0.8, ymax = 1) +
  annotate('richtext', x = 0.22, y = 1.01, label = "<span style='font-size:8pt'>A</span>", fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_initial_profiles$C), xmin = 0.3, xmax = 0.45, ymin = 0.8, ymax = 1) +
  annotate('richtext', x = 0.38, y = 1.01, label = "<span style='font-size:8pt'>C</span>", fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_initial_profiles$D), xmin = 0.46, xmax = 0.61, ymin = 0.8, ymax = 1) +
  annotate('richtext', x = 0.54, y = 1.01, label = "<span style='font-size:8pt'>D</span>", fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_initial_profiles$W), xmin = 0.66, xmax = 0.81, ymin = 0.8, ymax = 1) +
  annotate('richtext', x = 0.74, y = 1.01, label = "<span style='font-size:8pt'>W</span>", fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_initial_profiles$Y), xmin = 0.81, xmax = 0.96, ymin = 0.8, ymax = 1) +
  annotate('richtext', x = 0.89, y = 1.01, label = "<span style='font-size:8pt'>Y</span>", fill = NA, label.colour = NA) +
  annotate('point', x = c(0.625, 0.6425, 0.66), y = 0.9, shape = 20, size = 0.5) +
  
  # Connecting arrows
  annotate('segment', x = 0.225, xend = 0.225, y = 0.84, yend = 0.83) +
  annotate('curve', x = 0.225, xend = 0.26, y = 0.83, yend = 0.795, curvature = 0.35) +
  annotate('curve', x = 0.26, xend = 0.295, y = 0.795, yend = 0.755, curvature = -0.35) +
  annotate('segment', x = 0.295, xend = 0.295, y = 0.75, yend = 0.745,
           arrow.fill = 'black', arrow = arrow(type = 'closed', length = unit(0.02, 'npc'))) +
  annotate('curve', x = 0.225, xend = 0.19, y = 0.83, yend = 0.795, curvature = -0.35) +
  annotate('curve', x = 0.19, xend = 0.155, y = 0.795, yend = 0.755, curvature = 0.35) +
  annotate('segment', x = 0.155, xend = 0.155, y = 0.75, yend = 0.745,
           arrow.fill = 'black', arrow = arrow(type = 'closed', length = unit(0.02, 'npc'))) +
  
  # Permssive / Other A positions
  annotation_custom(ggplotGrob(p_a_profiles$AP), xmin = 0.075, xmax = 0.225, ymin = 0.445, ymax = 0.745) +
  annotate('richtext', x = 0.155, y = 0.7, label = "<span style='font-size:8pt'>AP</span>", fill = NA, label.colour = NA) +
  annotation_custom(ggplotGrob(p_a_profiles$Rest), xmin = 0.225, xmax = 0.375, ymin = 0.445, ymax = 0.745) +
  annotate('richtext', x = 0.305, y = 0.7, label = "<span style='font-size:8pt'>Rest</span>", fill = NA, label.colour = NA) +
  
  # Dendrograms and arrows
  annotation_custom(ggplotGrob(p_schematic_dend), xmin = 0, xmax = 0.45, ymin = 0.24, ymax = 0.5) +
  annotate('curve', x = 0.305, xend = 0.27, y = 0.525, yend = 0.5, curvature = -0.35) +
  annotate('curve', x = 0.27, xend = 0.235, y = 0.5, yend = 0.465, curvature = 0.35) +
  annotate('segment', x = 0.235, xend = 0.235, y = 0.457, yend = 0.452,
           arrow.fill = 'black', arrow = arrow(type = 'closed', length = unit(0.02, 'npc'))) +
  
  annotation_custom(ggplotGrob(p_schematic_dend_cuts), xmin = 0, xmax = 0.45, ymin = 0, ymax = 0.26) +
  annotate('segment', x = 0.235, xend = 0.235, y = 0.2875, yend = 0.2125, arrow.fill = 'black', arrow = arrow(type = 'closed', length = unit(0.02, 'npc'))) +
  
  # Text
  annotate('TextBox', x = 0.405, y = 0.825, hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>1.<br></span>") +
  annotate('TextBox', x = 0.455, y = 0.825, width = unit(0.65, 'npc'), hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>Process each AA independantly<br></span>") +
  annotate('TextBox', x = 0.405, y = 0.665, hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>2.<br></span>") +
  annotate('TextBox', x = 0.455, y = 0.675, width = unit(0.65, 'npc'), hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>Split permissive positions (|ER|&nbsp;<&nbsp;0.4)<br></span>") +
  annotate('TextBox', x = 0.405, y = 0.5, hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>3.<br></span>") +
  annotate('TextBox', x = 0.455, y = 0.5, width = unit(0.65, 'npc'), hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>Heirarchical clustering using cosine distance between PC2:20 profiles<br></span>") +
  annotate('TextBox', x = 0.405, y = 0.25, hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>4.<br></span>") +
  annotate('TextBox', x = 0.455, y = 0.25, width = unit(0.65, 'npc'), hjust = 0, vjust = 1, fill = NA, box.colour = NA, lineheight = 0.6,
           label = "<span style='font-size:7pt'>Determines subtypes with hybrid dynamic tree cutting<br></span>")

### Panel 2 - Cluster Sizes ###
subtype_size <- group_by(dms, cluster, wt) %>%
  summarise(size = n()) %>%
  group_by(wt) %>%
  mutate(n_aa = sum(size)) %>%
  ungroup() %>%
  mutate(freq = size / n_aa,
         cluster_num = str_sub(cluster, -1),
         wt = factor(wt, levels = rev(sort(unique(wt)))),
         cluster_num = factor(cluster_num, levels = c('O', 'P', 8:1)))

aa_labs <- levels(subtype_size$wt)
n_aa_labs <- structure(subtype_size$n_aa, names = as.character(subtype_size$wt))[aa_labs]

p_sizes <- ggplot(subtype_size, aes(y = as.integer(wt), x = freq, fill = cluster_num)) +
  geom_col() +
  scale_y_continuous(breaks = 1:length(aa_labs), labels = aa_labs,
                     sec.axis = sec_axis(~., breaks = 1:length(n_aa_labs), labels = n_aa_labs, name = ' Total Positions'), expand = expansion(0, 0)) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_fill_manual(values = CLUSTER_NUM_COLOURS) +
  guides(fill = guide_legend(title = 'Subtype', reverse = TRUE, nrow = 1, label.position = 'top', title.hjust = 0.5,
                             direction = 'horizontal', title.position = 'top')) +
  labs(x = 'Frequency') +
  theme(axis.title.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = 'top',
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(3, 'mm'),
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(0, 0, -2, 0, 'mm'),
        legend.box.background = element_blank(),
        legend.background = element_blank())

### Panel 3 - Correlation plot/dendrogram with labels ###
cor_subtype_sets <- list(`Small Aliphatic` = c('G1', 'G2', 'A5', 'G3', 'S1'),
                         `Not Proline` = c('I3', 'Y4', 'T2', 'E2', 'D3', 'L6', 'N2', 'R2', 'K3', 'A3', 'V5', 'S2', 'Q2', 'M2'),
                         `Positive` = c('Q4', 'K5', 'R3', 'K2', 'R5', 'K1', 'R1'),
                         Aromatic = c('Y3', 'Y1', 'H1', 'W1', 'F2'),
                         Aliphatic = c('C2', 'A2', 'G7', 'S3', 'V1', 'L3', 'I1', 'P4', 'L4'),
                         `Large Aliphatic` = c('T4', 'L5', 'L1', 'F1', 'Y2', 'V3', 'V2', 'M1', 'I2', 'L2'),
                         `Not Aromatic` = c('T6', 'P2', 'T1', 'A1', 'P3', 'A4', 'T3', 'R4'),
                         `Negative` = c('T5', 'Q3', 'N1', 'Q1', 'E3', 'L7', 'G4', 'D2', 'D1', 'E1'))

cor_set_colours <- c(`Small Aliphatic` = '#ff7f00', Aliphatic = '#ffff33', `Large Aliphatic` = '#a65628',
                     `Not Proline` = '#4daf4a', Permissive = '#999999',
                     `Positive` = '#984ea3', `Negative` = '#f781bf',
                     Aromatic = '#377eb8', `Not Aromatic` = '#e41a1c') # Manual assignment of colourbrewer2 Set1 colours

classify_cluster <- function(x){
  out <- rep(NA, length(x))
  out[str_detect(x, CLUSTER_PERMISSIVE_RE)] <- 'Permissive'
  for (set in names(cor_subtype_sets)){
    out[x %in% cor_subtype_sets[[set]]] <- set
  }
  return(out)
}

# Correlation heatmap/dendrogram
cors <- filter(dms, !str_detect(cluster, CLUSTER_OUTLIER_RE), !str_detect(cluster, CLUSTER_PERMISSIVE_RE)) %>%
  cluster_profile_correlation(A:Y)

hc <- select(cors, cluster1, cluster2, cor) %>%
  pivot_wider(id_cols = cluster1, names_from = cluster2, values_from = cor) %>%
  tibble_to_matrix(-cluster1, row_names = 'cluster1') %>%
  subtract(1, .) %>%
  as.dist() %>%
  hclust()
dend_data <- dendro_data(hc)
branches <- dend_data$segments
leaves <- dend_data$labels %>%
  mutate(type = classify_cluster(label))

cors <- mutate(cors, cluster1 = factor(cluster1, levels = leaves$label), cluster2 = factor(cluster2, levels = leaves$label))

annotate_set_rect <- function(set){
  inds <- which(levels(cors$cluster1) %in% cor_subtype_sets[[set]])
  mn <- min(inds) - 0.5
  mx <- max(inds) + 0.5
  x = c(mn, mx, mx, mn)
  y = c(mn, mn, mx, mx)
  annotate('polygon', x = x, y = y, colour = cor_set_colours[set], fill = NA)
}

p_cor_heatmap <- ggplot(cors, aes(x=as.integer(cluster1), y=as.integer(cluster2), fill=cor)) +
  geom_tile() +
  geom_segment(data = branches, aes(x=-5*y-1, y=x, xend=-5*yend-1, yend=xend), inherit.aes = FALSE) +
  geom_point(data = leaves, aes(x=-5*y-1, y=x, colour = type), shape=20, size = 0.8, inherit.aes = FALSE) +
  annotate_set_rect("Small Aliphatic") +
  annotate_set_rect("Not Proline") +
  annotate_set_rect("Positive") +
  annotate_set_rect("Aromatic") +
  annotate_set_rect("Aliphatic") +
  annotate_set_rect("Large Aliphatic") +
  annotate_set_rect("Not Aromatic") +
  annotate_set_rect("Negative") +
  scale_fill_distiller(type = ER_COR_COLOURS$type, palette = ER_COR_COLOURS$palette, direction = ER_COR_COLOURS$direction,
                       limits = c(-1, 1)) +
  scale_colour_manual(values = cor_set_colours, na.value='grey') +
  coord_fixed() +
  guides(fill = guide_colourbar(title = 'Pearson\nCorrelation', barwidth = unit(3, 'mm'), barheight = unit(15, 'mm')), colour = FALSE) +
  theme(text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'mm'),
        legend.background = element_blank())
p_cor_heatmap_legend <- as_ggplot(get_legend(p_cor_heatmap))
p_cor_heatmap <- p_cor_heatmap + guides(fill = FALSE)

# Profiles of key groups
plot_profile_block <- function(set, guide=FALSE){
  p <- filter(full_characterisation$profiles, cluster %in% set) %>%
    mutate(cluster = factor(cluster, levels = levels(leaves$label)[levels(leaves$label) %in% set]),
           mut = factor(mut)) %>%
    ggplot(aes(x = mut, y = cluster, fill = er)) +
    geom_tile(colour = 'grey', size = 0.05) +
    coord_fixed() +
    scale_fill_distiller(type = ER_PROFILE_COLOURS$type, palette = ER_PROFILE_COLOURS$palette, direction = ER_PROFILE_COLOURS$direction,
                         limits = c(-1, 1)) +
    theme(text = element_text(size = 6),
          plot.title = element_text(margin = margin(0, 0, 0, 0, 'mm')),
          axis.ticks = element_blank(),
          axis.text.y = element_markdown(margin = margin(0, 0, 0, 0, 'mm')),
          axis.text.x = element_markdown(margin = margin(0, 0, 0, 0, 'mm')),
          panel.background = element_blank(),
          axis.title = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.background = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'mm'))
  if (guide){
    p <- p + guides(fill = guide_colourbar(title = 'Mean ER', barwidth = unit(3, 'mm'), barheight = unit(15, 'mm')))
  } else {
    p <- p + guides(fill = FALSE)
  }
  p
}

prof_title <- function(x){
  str_c("<span style='font-size:12pt; color:", cor_set_colours[x],"'>&#9679;</span> ", x)
}

p_cor_set_profiles <- map(names(cor_subtype_sets), ~plot_profile_block(cor_subtype_sets[[.]]) + 
                            labs(title = prof_title(.)) + 
                            theme(plot.title = element_markdown(size = 7))) %>%
  set_names(names(cor_subtype_sets))
p_cor_profiles_legend <- plot_profile_block(cor_subtype_sets$small_aliphatic, guide = TRUE) %>% 
  get_legend() %>% 
  as_ggplot()

## Main plot
p_heatmap <- ggplot() +
  geom_blank() +
  lims(x = c(0, 1), y = c(0, 1)) + 
  coord_fixed(clip = 'off') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank()) +

  annotation_custom(ggplotGrob(p_cor_heatmap), xmin = 0, xmax = 0.75, ymin = 0.13, ymax = 0.87) +
  annotation_custom(ggplotGrob(p_cor_heatmap_legend), xmin = -0.075, xmax = 0.075, ymin = 0.5, ymax = 0.75) +
  annotation_custom(ggplotGrob(p_cor_profiles_legend), xmin = -0.1, xmax = 0.075, ymin = 0.25, ymax = 0.5) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Positive`), xmin = -0.05, xmax = 0.3, ymin = 0.79, ymax = 1.03) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Small Aliphatic`), xmin = 0.31, xmax = 0.71, ymin = 0.8, ymax = 1) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Not Proline`), xmin = 0.72, xmax = 1.07, ymin = 0.65, ymax = 1) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$Aromatic), xmin = 0.72, xmax = 1.07, ymin = 0.475, ymax = 0.675) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$Aliphatic), xmin = 0.72, xmax = 1.07, ymin = 0.25, ymax = 0.47) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Large Aliphatic`), xmin = 0.72, xmax = 1.07, ymin = 0, ymax = 0.25) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Not Aromatic`), xmin = 0.3, xmax = 0.7, ymin = 0, ymax = 0.2) +
  annotation_custom(ggplotGrob(p_cor_set_profiles$`Negative`), xmin = -0.1, xmax = 0.33, ymin = 0, ymax = 0.2)

### Panel 4 - Clusters mapped to UMAP ###
cluster_dms <- mutate(dms, cluster_type = classify_cluster(cluster)) %>%
  select(gene, position, wt, umap1, umap2, cluster, cluster_type) %>%
  drop_na(cluster_type)

p_umap <- ggplot() +
  geom_point(data = dms, mapping = aes(x = umap1, y = umap2), colour = 'grey90', shape = 20, size = 1) +
  geom_point(data = cluster_dms, mapping = aes(x = umap1, y = umap2, colour = cluster_type), shape = 20, size = 1) +
  scale_color_manual(values = cor_set_colours) +
  labs(x = 'UMAP1', y = 'UMAP2') + 
  guides(colour = guide_legend(title = '', nrow = 3, override.aes = list(size = 2))) +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.key.height = unit(4, 'mm'),
        legend.key.width = unit(2, 'mm'),
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(0, 0, -5, -10, 'mm'),
        legend.box.background = element_blank(),
        legend.background = element_blank())

### Panel 5 - Subtype frequencies ###
freq_summary <- group_by(full_characterisation$summary, aa) %>%
  mutate(freq = n / sum(n)) %>%
  summarise(Permissive = freq[which(cluster == str_c(aa, 'P'))],
            `Not Proline` = max(freq[which(cluster %in% cor_subtype_sets$`Not Proline`)], 0),
            `Most Selective` = freq[which(cluster %in% most_selective_subtypes$cluster)],
            Other = 1 - (Permissive + `Not Proline` + `Most Selective`)) %>%
  pivot_longer(-aa, names_to = 'type', values_to = 'freq') %>%
  mutate(type = factor(type, levels = c('Permissive', 'Not Proline', 'Other', 'Most Selective'))) %>%
  left_join(select(most_selective_subtypes, aa, mean_er), by = 'aa')

p_subtype_freqs <- ggplot(freq_summary) +
  geom_col(aes(x = freq, y = reorder(aa, desc(aa)), fill = type)) +
  scale_fill_brewer(type = 'qual', palette = 'Paired') + 
  guides(fill = guide_legend(title = '', direction = 'horizontal', reverse = TRUE, label.position = 'top')) +
  labs(x = 'Frequency') +
  theme(panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'top',
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(3, 'mm'),
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(0, 0, -2, 0, 'mm'),
        legend.box.background = element_blank(),
        legend.background = element_blank())

### Panel 6 - Profile of most selective subtype ###
most_selective_profiles <- filter(full_characterisation$profiles, cluster %in% most_selective_subtypes$cluster) %>%
  mutate(wt = str_sub(cluster, end = 1))

p_selective_profiles <- ggplot(most_selective_profiles, aes(y = reorder(wt, desc(wt)), x = mut, fill = er)) +
  geom_tile(colour = 'grey') +
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = -1, limits = c(-0.8, max(most_selective_profiles$er))) +
  labs(x = 'Substitution', y = '') +
  guides(fill = guide_colourbar(title = 'ER', reverse = TRUE, direction = 'horizontal', title.vjust = 0,
                                barheight = unit(2.5, 'mm'), label.position = 'top')) +
  theme(axis.ticks = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.1, 'lines'),
        panel.grid.major.y = element_blank(),
        legend.position = 'top',
        legend.margin = margin(0, 0, 0, 0, 'mm'),
        legend.box.margin = margin(0, 0, -2, -5, 'mm'),
        legend.box.background = element_blank(),
        legend.background = element_blank())

### Assemble figure ###
p1 <- p_schematic + labs(tag = 'A') + theme(text = element_text(size = 8))
p2 <- p_sizes + labs(tag = 'B') + theme(text = element_text(size = 7.5))
p3 <- p_heatmap + labs(tag = 'C') + theme(text = element_text(size = 8.5))
p4 <- p_umap + labs(tag = 'D') + theme(text = element_text(size = 8))
p5 <- p_subtype_freqs + labs(tag = 'E') + theme(text = element_text(size = 8))
p6 <- p_selective_profiles + labs(tag = 'F') + theme(text = element_text(size = 8.5))

figure3 <- multi_panel_figure(width = 183, height = 183, columns = 9, rows = 3,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:3) %>%
  fill_panel(p2, row = 2, column = 1:3) %>%
  fill_panel(p3, row = 1:2, column = 4:9) %>%
  fill_panel(p4, row = 3, column = 1:3) %>%
  fill_panel(p5, row = 3, column = 4:6) %>%
  fill_panel(p6, row = 3, column = 7:9)

ggsave('figures/4_figures/figure3.pdf', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm', device = cairo_pdf)
ggsave('figures/4_figures/figure3.png', figure3, width = figure_width(figure3), height = figure_height(figure3), units = 'mm')
