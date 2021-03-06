#!/usr/bin/env Rscript
# Check the validity of the normalisation approach
source('src/config.R')

dms <- read_tsv('data/long_combined_mutational_scans.tsv')
plots <- list()

### Check normalisation 10% with an without nonsense ###
quantiles <- filter(dms, !is.na(score)) %>%
  group_by(study) %>%
  summarise(with = quantile(transformed_score, 0.1, na.rm = TRUE),
            without = quantile(transformed_score[!mut == "*"], 0.1, na.rm = TRUE),
            nonsense = any(mut == "*"))

plots$norm_quantiles <- ggplot(filter(quantiles, nonsense), aes(x = with, y = without)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "With Nonsense", y = "Without Nonsense")

### Check distribution of bottom 10% scores predictor metrics ###
dms_scores <- select(dms, study, score, transformed_score, sift, log10_sift, total_energy) %>%
  left_join(quantiles, by = "study") %>%
  mutate(with10 = transformed_score < with,
         without10 = transformed_score < without,
         neutral = abs(score) < 0.3)

score_groups <- filter(dms_scores, nonsense, with10 | without10 | neutral) %>%
  pivot_longer(c(log10_sift, total_energy), names_to = "metric", values_to = "value") %>%
  pivot_longer(c(with10, without10, neutral), names_to = "class", values_to = "member") %>%
  filter(member, !is.na(value))

plots$score_group_density <- ggplot(score_groups, aes(x = value, y = ..scaled.., colour = class)) +
  stat_density(geom = "line", position = "identity") +
  facet_wrap(~metric, nrow = 1, scales = "free_x")

plots$score_group_box <- ggplot(score_groups, aes(x = class, y = value, fill = class)) +
  geom_boxplot() +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  stat_compare_means(comparisons = list(c("with10", "without10"), c("with10", "neutral"), c("without10", "neutral")))

### SIFT4G Scores across studies

plots$per_study_sift_dist <- ggplot(filter(dms_scores, with10), aes(x = study, y = log10_sift)) +
  geom_boxplot(fill = "cornflowerblue") +
  geom_hline(yintercept = log10(0.05)) +
  coord_flip() +
  labs(y = expression("log"[10]~"SIFT4G"), x = "") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"))

plots$per_study_foldx_dist <- ggplot(filter(dms_scores, with10), aes(x = study, y = total_energy)) +
  geom_boxplot(fill = "cornflowerblue") +
  geom_hline(yintercept = c(-1, 1)) +
  coord_flip() +
  labs(y = expression(Delta*Delta*"G"), x = "") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"))

### Save plots ###
save_plotlist(plots, "figures/0_data/normalisation", overwrite = "all")
