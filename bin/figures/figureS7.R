#!/usr/bin/env Rscript
# Produce figure S7 (Frequency of rarer subtypes)
source('src/config.R')
source('src/subtype_characterisation.R')

dms <- full_join(read_tsv('data/subtypes/final_subtypes.tsv'),
                 read_tsv('data/combined_mutational_scans.tsv'),
                 by = c('study', 'gene', 'position', 'wt')) %>%
  arrange(study, position)

full_characterisation <- full_cluster_characterisation(dms)
n_clusters <- nrow(full_characterisation$summary)

### Panel 1 - Frequencies ###
sub_freqs <- group_by(full_characterisation$summary, aa) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  filter(!str_detect(cluster, CLUSTER_PERMISSIVE_RE), !str_detect(cluster, CLUSTER_OUTLIER_RE)) %>% 
  mutate(num = as.integer(str_sub(cluster, 2)))

model <- lm(log(freq) ~ num + 0, data = sub_freqs)
coefs <- summary(model)$coefficients['num',]
x <- c(seq(1, 8, 0.1), seq(8, 1, -0.1))
y <- c(exp((coefs['Estimate']+coefs['Std. Error'])*seq(1, 8, 0.1)),
       exp((coefs['Estimate']-coefs['Std. Error'])*seq(8, 1, -0.1)))

p_freq <- ggplot(sub_freqs, aes(x = num, y = freq)) + 
  geom_point(aes(colour = aa)) +
  annotate('polygon', x = x, y = y, fill = 'lightgrey', alpha = 0.5, colour = NA) +
  stat_function(fun = ~exp(coefs['Estimate']*.)) +
  scale_colour_manual(values = AA_COLOURS) +
  scale_x_continuous(breaks = 1:8) +
   guides(colour = guide_legend(title = '', direction = 'horizontal', nrow = 2, byrow = TRUE)) +
   labs(x = 'Subtype', y = 'Frequency') +
   theme(legend.position = 'bottom')

### Panel 2 - Number of subtypes ###
# Generated by bin/analysis/2_subtypes/estimate_number.R
counts <- read_tsv("data/subtypes/rep_counts.tsv")

count_means <- group_by(counts, rep, positions) %>% 
  summarise(nclusters = sum(nclusters), .groups = "drop") %>%
  bind_rows(`All Amino Acids` = .,
            A = filter(counts, wt == "A"),
            C = filter(counts, wt == "C"),
            H = filter(counts, wt == "H"),
            W = filter(counts, wt == "W"),
            .id = "type") %>% 
  select(type, rep, positions, nclusters) %>%
  group_by(type, positions) %>%
  summarise(mean = mean(nclusters),
            sd = sd(nclusters),
            median = median(nclusters),
            .groups = "drop")

p_count_all <- ggplot(filter(count_means, type == "All Amino Acids"), aes(x = positions)) +
  geom_point(aes(y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000, 5000, 6000, 6357)) +
  labs(x = "Positions Clustered", y = "Functional Subtypes") +
  lims(y = c(0, 100))

p_count_aas <- ggplot(filter(count_means, type != "All Amino Acids"), aes(x = positions, colour = type)) +
  geom_point(aes(y = mean), show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), show.legend = FALSE) +
  scale_colour_manual(name = "", values = AA_COLOURS) +
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000, 5000, 6000, 6357)) +
  labs(x = "Positions Clustered", y = "Functional Subtypes") +
  theme(legend.position = "top")

### Assemble Figure ###
size <- theme(text = element_text(size = 8))
p1 <- p_freq + labs(tag = 'A') + size
p2 <- p_count_all + labs(tag = 'B') + size
p3 <- p_count_aas + labs(tag = 'C') + size

figure <- multi_panel_figure(width = 183, height = 183, columns = 2, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:2) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 2, column = 2)

ggsave('figures/4_figures/figureS7.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS7.png', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS7.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS7.eps', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device=cairo_ps, fallback_resolution = 600)

