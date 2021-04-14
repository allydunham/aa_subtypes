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

plots$norm_quantiles <- ggplot(quantiles[quantiles$nonsense], aes(x = with, y = without)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "With Nonsense", y = "Without Nonsense")

### Check distribution of bottom 10% scores predictor metrics ###
dms_scores <- select(dms, study, transformed_score, log10_sift, total_energy) %>%
  left_join(quantiles, by = "study") %>%
  mutate(with10 = transformed_score < with,
         without10 = transformed_score < without)

bot10 <- filter(dms_scores, nonsense, with10 | without10) %>%
  pivot_longer(c(log10_sift, total_energy), names_to = "metric", values_to = "value") %>%
  pivot_longer(ends_with("10"), names_to = "class", values_to = "member") %>%
  filter(member, !is.na(value))

plots$bottom_10_scores <- ggplot(bot10, aes(x = value, y = ..scaled.., colour = class)) +
  stat_density(geom = "line", position = "identity") +
  facet_wrap(~metric, nrow = 1, scales = "free_x")
