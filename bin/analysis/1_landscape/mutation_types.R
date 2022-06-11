#!/usr/bin/env Rscript
# Analyse mutation types looking for patterns
source('src/config.R')
library(caret)

dms <- read_tsv('data/combined_mutational_scans.tsv') %>%
  mutate(uniprot_id = unname(UNIPROT_IDS[gene]))
  
dms_correlation <- tibble_correlation(dms, A:Y)

p_cor <- ggplot(dms_correlation, aes(x = cat1, y = cat2, fill = cor)) +
  geom_raster() +
  scale_fill_distiller(name = expression(rho), palette = "YlOrRd", limits = c(0, 1), direction = 1) +
  coord_equal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(hjust = 0.15),
        legend.key.height = unit(30, "pt"))

set.seed(1995)
trainIndex <- createDataPartition(dms$A, p = 0.95, times = 1, list = FALSE)
training <- dms[trainIndex,]
testing  <- dms[-trainIndex,]

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

# Linear model
models <- list(
  A = train(A ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  E = train(E ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  H = train(H ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  L = train(L ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  M = train(M ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  N = train(N ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  Q = train(Q ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  R = train(R ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  S = train(S ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  V = train(V ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  W = train(W ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl),
  Y = train(Y ~ `F` + I + C + `T` + P + K + G + D, data = dms, method = "lm", trControl = fitControl)
)

model_performance <- map_dfr(models, ~.$results, .id = "aa") %>%
  select(-intercept) %>%
  pivot_longer(-aa, names_to = "metric", values_to = "value") %>%
  mutate(metric = ifelse(str_ends(metric, "SD"), metric, str_c(metric, "value"))) %>%
  tidyr::extract(metric, c("metric", "type"), "([A-Za-z]*)(SD|value)") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  rename_with(str_to_lower)

p_models <- ggplot(model_performance, aes(x = aa, y = value, ymin = value - sd, ymax = value + sd, fill = metric)) +
  facet_wrap(~metric, scales = "free_y", ncol = 1, strip.position = "left") +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(width = 0.3) +
  scale_fill_brewer(palette = "Dark2") +
  theme(strip.placement = "outside",
        axis.title = element_blank())

# Classifier model
classifier_train <- mutate(dms, across(A:Y, .fns = list(class = ~as.factor(. < -0.5)))) %>%
  select(A:Y, A_class:Y_class)

classifier_models <- list(
  A = train(A_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  E = train(E_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  H = train(H_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  L = train(L_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  M = train(M_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  N = train(N_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  Q = train(Q_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  R = train(R_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  S = train(S_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  V = train(V_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  W = train(W_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl),
  Y = train(Y_class ~ `F` + I + C + `T` + P + K + G + D, data = classifier_train, method = "glm", family = "binomial", trControl = fitControl)
)

classifier_performance <- map_dfr(classifier_models, ~.$results, .id = "aa") %>%
  select(-parameter) %>%
  pivot_longer(-aa, names_to = "metric", values_to = "value") %>%
  mutate(metric = ifelse(str_ends(metric, "SD"), metric, str_c(metric, "value"))) %>%
  tidyr::extract(metric, c("metric", "type"), "([A-Za-z]*)(SD|value)") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  rename_with(str_to_lower)

p_classifier <- ggplot(classifier_performance, aes(x = aa, y = value, ymin = value - sd, ymax = value + sd, fill = metric)) +
  facet_wrap(~metric, scales = "free_y", ncol = 1, strip.position = "left") +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(width = 0.3) +
  scale_fill_brewer(palette = "Dark2") +
  theme(strip.placement = "outside",
        axis.title = element_blank())
