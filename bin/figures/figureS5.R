#!/usr/bin/env Rscript
# Produce figure S5 (Multiple variants)
source('src/config.R')
source('src/study_standardising.R')

# Araya 2012
read_araya <- function(){
  meta <- read_yaml('data/studies/araya_2012_yap1/araya_2012_yap1.yaml')
  dm_data <- read_csv('data/studies/araya_2012_yap1/raw/urn_mavedb_00000002-a-2_scores.csv', skip = 4) %>%
    select(hgvs_pro, raw_score = score) %>%
    mutate(hgvs_pro = if_else(str_ends(hgvs_pro, ']'), str_sub(hgvs_pro, start = 4, end = -2), str_sub(hgvs_pro, start = 3)),
           n_mut = str_count(hgvs_pro, ';') + 1) %>%
    separate(hgvs_pro, str_c('mut', 1:max(.$n_mut)), sep = ';', fill = 'right') %>%
    pivot_longer(cols = starts_with('mut'), values_to = 'mut') %>%
    drop_na(mut) %>%
    select(-name) %>%
    tidyr::extract(mut, into = c('wt', 'position', 'mut'), "([A-Za-z]{3})([0-9]+)([A-Za-z]{3})", convert = TRUE) %>%
    mutate(wt = AA_THREE_2_ONE[wt], mut = AA_THREE_2_ONE[mut], position = position + 169) %>%
    mutate(transformed_score = raw_score) %>%
    group_by(position, wt, mut)
  
  singles <- summarise(dm_data, single_score := ifelse(any(n_mut == 1), mean(transformed_score[n_mut == 1], na.rm = TRUE), NA)) %>%
    ungroup()
  single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)
  
  lapply(2:9, function(x){
    summarise(dm_data, !!str_c('score_', x) := ifelse(any(n_mut <= x), mean(transformed_score[n_mut <= x], na.rm = TRUE), NA)) %>%
      ungroup() %>%
      select(!!str_c('score_', x))
  }) %>%
    bind_cols(singles, .) %>%
    pivot_longer(starts_with('score_'), names_to = 'n_mut', values_to = 'multi_score', names_prefix = 'score_', names_ptypes = list(n_mut=integer())) %>%
    group_by(n_mut) %>%
    mutate(frac = sum(!is.na(multi_score))/length(multi_score)) %>%
    ungroup() %>%
    mutate(n_mut = str_c(n_mut, ' (', signif(frac, digits = 4)*100, '%)'),
           study = str_c('Araya et al. 2012 (YAP1)\n', signif(single_frac, digits = 4)*100, '% measures as single variants'))
}

# Dorrity 2018
read_dorrity <- function(){
  meta <- read_yaml('data/studies/dorrity_2018_ste12/dorrity_2018_ste12.yaml')
  mating_data <- read_xlsx('data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd01.xlsx') %>%
    mutate(mut = sapply(seqID, process_split_seqid, USE.NAMES = FALSE),
           n_mut = str_count(mut, ',') + 1,
           mating_avg = (mating_30C_rep1 + mating_30C_rep2 + mating_30C_rep3)/3) %>%
    select(mut, n_mut, mating_avg, starts_with('mating_30C_'))
  
  invasion_data <- read_xlsx('data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd02.xlsx') %>%
    mutate(mut = sapply(seqID, process_split_seqid, USE.NAMES = FALSE),
           n_mut = str_count(mut, ',') + 1,
           invasion_avg = (invasion_30C_rep1 + invasion_30C_rep2 + invasion_30C_rep3)/3) %>%
    select(mut, n_mut, invasion_avg, starts_with('invasion_30C_'))
  
  dm_data <- full_join(mating_data, invasion_data, by = c('mut', 'n_mut')) %>%
    mutate(raw_score = pmin(mating_avg, invasion_avg, na.rm = TRUE)) %>%
    select(mut, n_mut, raw_score) %>%
    separate(mut, into = str_c('mut', 1:max(.$n_mut)), sep=',', fill = 'right') %>%
    pivot_longer(starts_with('mut'), names_to = 'tmp', values_to = 'mut') %>%
    select(-tmp) %>%
    drop_na(mut) %>%
    tidyr::extract(mut, into = c('position', 'mut'), '([0-9]*)([A-Z*])', convert = TRUE) %>%
    mutate(position = position + 140,
           wt = str_split(meta$seq, '')[[1]][position])
  
  singles <- group_by(dm_data, position, wt, mut) %>%
    summarise(single_score := ifelse(any(n_mut == 1), mean(raw_score[n_mut == 1], na.rm = TRUE), NA)) %>%
    ungroup()
  single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)
  
  lapply(2:max(dm_data$n_mut), function(x){
    group_by(dm_data, position, wt, mut) %>%
      summarise(!!str_c('score_', x) := ifelse(any(n_mut <= x), mean(raw_score[n_mut <= x], na.rm = TRUE), NA)) %>%
      ungroup() %>%
      select(!!str_c('score_', x))
  }) %>%
    bind_cols(singles, .) %>%
    pivot_longer(starts_with('score_'), names_to = 'n_mut', values_to = 'multi_score', names_prefix = 'score_', names_ptypes = list(n_mut=integer())) %>%
    group_by(n_mut) %>%
    mutate(frac = sum(!is.na(multi_score))/length(multi_score)) %>%
    ungroup() %>%
    mutate(n_mut = str_c(n_mut, ' (', signif(frac, digits = 4)*100, '%)'),
           study = str_c('Dorrity et al. 2018 (STE12)\n', signif(single_frac, digits = 4)*100, '% measures as single variants'))
}

# Starita 2013
read_starita <- function(){
  meta <- read_yaml('data/studies/starita_2013_ube4b/starita_2013_ube4b.yaml')
  dm_data <- read_xlsx('data/studies/starita_2013_ube4b/raw/sd01.xlsx', na = c('NA', '')) %>%
    filter(!seqID == 'NA-NA') %>% # Filter WT
    rename(raw_score = log2_ratio) %>%
    separate(seqID, into = c('position', 'mut'), sep='-') %>%
    select(-nscor_log2_ratio) %>%
    mutate(n_mut = sapply(position, function(x){str_count(x, ',') + 1})) %>%
    separate(mut, str_c('mut', 1:max(.$n_mut)), sep = ',', fill = 'right') %>%
    separate(position, str_c('position', 1:max(.$n_mut)), sep = ',', fill = 'right') %>%
    pivot_longer(starts_with('position'), names_to = 'pos_num', names_prefix = 'position', values_to = 'position') %>%
    drop_na(position) %>%
    pivot_longer(starts_with('mut'), names_to = 'mut_num', names_prefix = 'mut', values_to = 'mut') %>%
    drop_na(mut) %>%
    filter(pos_num == mut_num) %>%
    select(-pos_num, -mut_num) %>%
    group_by(position, mut)
  
  singles <- summarise(dm_data, single_score := ifelse(any(n_mut == 1), mean(raw_score[n_mut == 1], na.rm = TRUE), NA)) %>%
    ungroup()
  single_frac <- sum(!is.na(singles$single_score))/length(singles$single_score)
  
  lapply(2:9, function(x){
    summarise(dm_data, !!str_c('score_', x) := ifelse(any(n_mut <= x), mean(raw_score[n_mut <= x], na.rm = TRUE), NA)) %>%
      ungroup() %>%
      select(!!str_c('score_', x))
  }) %>%
    bind_cols(singles, .) %>%
    pivot_longer(starts_with('score_'), names_to = 'n_mut', values_to = 'multi_score', names_prefix = 'score_', names_ptypes = list(n_mut=integer())) %>%
    group_by(n_mut) %>%
    mutate(frac = sum(!is.na(multi_score))/length(multi_score)) %>%
    ungroup() %>%
    mutate(n_mut = str_c(n_mut, ' (', signif(frac, digits = 4)*100, '%)'),
           position = as.integer(position),
           study = str_c('Starita et al. 2013 (UBE4B)\n', signif(single_frac, digits = 4)*100, '% single variant coverage'))
}

# Combine and plot
mut_count_data <- list(read_araya(), read_dorrity(), read_starita()) %>%
  bind_rows() %>%
  mutate(n_mut = factor(n_mut,
                        levels = unique(n_mut)[unique(n_mut) %>% 
                                                 str_match('([0-9]*) \\([0-9\\.]*%\\)') %>% 
                                                 extract(,2) %>% 
                                                 as.integer %>% 
                                                 order()]))


panels <- group_by(mut_count_data, study) %>% 
  group_map(~ggplot(., aes(x = single_score, y = multi_score)) +
              facet_wrap(~n_mut, ncol = 5) +
              coord_cartesian(clip = 'off') +
              geom_point(colour = 'cornflowerblue', shape = 20) + 
              geom_abline(slope = 1, linetype='dashed') +
              labs(x = 'Score (Single Variants)',
                   y = 'Score (Mean Over Multiple Variants)',
                   title = .y$study))

figure <- multi_panel_figure(width = 183, height = 270, columns = 1, rows = 3, unit = 'mm') %>%
  fill_panel(panels[[1]], row = 1, column = 1) %>%
  fill_panel(panels[[2]], row = 2, column = 1) %>%
  fill_panel(panels[[3]], row = 3, column = 1)
ggsave('figures/4_figures/figureS5.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/4_figures/figureS5.png', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
