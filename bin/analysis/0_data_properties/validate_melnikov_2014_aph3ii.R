#!/usr/bin/env Rscript
# Validate normalisation methodology for Melnikov et al. 2014 (APH(3')-II)
source('src/config.R')
source('src/study_standardising.R')

dir.create('figures/0_data_properties/per_study/melnikov_2014_aph3ii')

#### Import data ####
count_files <- grep('\\.aacounts\\.txt', dir('data/studies/melnikov_2014_aph3ii/raw/'), value = TRUE)
count_files <- count_files[!grepl('(S[12]\\_Ami|S3\\_Kan)', count_files)] # Discard bad runs - see 00README.txt from Melnikov et al. data

counts <- sapply(count_files, read_melnikov_table, simplify = FALSE) %>%
  set_names(gsub('(KKA2\\_|\\.aacounts\\.txt)', '', names(.))) # Set names to drug
bkg_counts <- counts[c('Bkg1', 'Bkg2')]
counts <- counts[which(!names(counts) %in% c('Bkg1', 'Bkg2'))]

dm_data <- mapply(melnikov_fitness, counts, names(counts), MoreArgs = list(bkg=bkg_counts), SIMPLIFY = FALSE) %>%
  bind_rows(.id = 'experiment') %>%
  separate(experiment, c('round', 'drug', 'library'), sep='_') %>%
  # Round and library contain the same information (plus round notes which needed a re-test at different MIC, which we alread accounted for) -> discard round
  select(position, wt, mut, score, drug, library) %>%
  mutate(rel_conc = 1/as.numeric(str_sub(drug, -1)), drug = str_sub(drug, 1, -3)) %>%
  pivot_wider(id_cols = c('position', 'wt', 'mut', 'drug', 'rel_conc'), names_from = library, values_from = 'score') %>%
  mutate(diff = abs(L1 - L2))
########

#### Filter libraries ####
# Correlation between the two libraries
p_lib_test <- ggplot(dm_data, aes(x=L1, y=L2, colour=drug)) +
  facet_grid(rows = vars(drug), cols = vars(rel_conc)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope = 1) +
  guides(colour=FALSE) +
  labs(title = "Relationship between duplicate libraries in Melnikov et al. 2014 (APH(3')-II)",
       subtitle = "Faceted by drug and drug concentration relative to MIC")
ggsave('figures/0_data_properties/per_study/melnikov_2014_aph3ii/initial_library_correlation.pdf', p_lib_test, units = 'cm', width = 18, height = 18)

# Filter conditions where libraries don't agree and observations where libraries differ a lot
# Then take average over libraries
dm_data <- filter(dm_data, !(drug == 'Ami' & rel_conc == 0.25), !(drug %in% c('G418', 'Ami', 'Kan') & rel_conc == 0.125)) %>%
  filter(diff < sd(diff, na.rm = TRUE) * 3) %>%
  mutate(score = (L1 + L2)/2) %>%
  drop_na(score)

# Filtered Correlation between the two libraries
p_lib_filtered <- ggplot(dm_data, aes(x=L1, y=L2, colour=drug)) +
  facet_wrap(~drug, nrow = 2) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope = 1) +
  guides(colour=FALSE) +
  labs(title = "Relationship between duplicate libraries in Melnikov et al. 2014 (APH(3')-II)",
       subtitle = "After bad filtering bad libraries and observations where L1 - L2 > 3sd")
ggsave('figures/0_data_properties/per_study/melnikov_2014_aph3ii/filtered_library_correlation.pdf', p_lib_filtered, units = 'cm', width = 18, height = 18)
########

#### Drug and Concentration Relationships ####
dm_data_rc <- select(dm_data, position, wt, mut, score, drug, rel_conc) %>%
  pivot_wider(names_from = rel_conc, values_from = score)
dm_data_rc <- bind_rows(`1 - 0.5`=select(dm_data_rc, drug, position, wt, mut, c1=`1`, c2=`0.5`),
                        `1 - 0.25`=select(dm_data_rc, drug, position, wt, mut, c1=`1`, c2=`0.25`),
                        `1 - 0.125`=select(dm_data_rc, drug, position, wt, mut, c1=`1`, c2=`0.125`),
                        `0.5 - 0.25`=select(dm_data_rc, drug, position, wt, mut, c1=`0.5`, c2=`0.25`),
                        `0.5 - 0.125`=select(dm_data_rc, drug, position, wt, mut, c1=`0.5`, c2=`0.125`),
                        `0.25 - 0.125`=select(dm_data_rc, drug, position, wt, mut, c1=`0.25`, c2=`0.125`),
                        .id = 'rel_conc')

p_rel_conc <- ggplot(dm_data_rc, aes(x=c1, y=c2, colour=drug)) +
  facet_grid(rows = vars(drug), cols = vars(rel_conc)) +
  geom_point() +
  geom_smooth(method='lm', colour = 'black') +
  geom_abline(slope = 1, linetype='dotted') +
  guides(colour=FALSE) +
  labs(title = "Relationship between drug concentrations in Melnikov et al. 2014 (APH(3')-II)",
       subtitle = "Concentrations relative to MIC of drug",
       x = 'log2(ER), Conc. 1', y = 'log2(ER), Conc. 2')
ggsave('figures/0_data_properties/per_study/melnikov_2014_aph3ii/rel_conc_correlation.pdf', p_rel_conc, units = 'cm', width = 25, height = 25)

# Behaviour changes a lot depending on concentration of drug - Use one?
p_drug_conc <- ggplot(dm_data, aes(x = score, y = ..density.., fill = drug)) +
  facet_grid(rows = vars(drug), cols = vars(rel_conc)) +
  geom_histogram() +
  geom_vline(xintercept = 0, linetype='dashed') +
  guides(fill=FALSE) +
  labs(title = "Distribution of log(ER) for each Drug-Conc. Combination",
       subtitle = "from Melnikov et al. 2014 (APH(3')-II)",
       x = 'log2(ER)', y = 'Frequency')
ggsave('figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_conc_distribution.pdf', p_drug_conc, units = 'cm', width = 25, height = 25)

# Choose best distributions:
dm_data <- filter(dm_data,
                  (drug == 'Ami' & rel_conc == 0.5) |
                    (drug == 'G418' & rel_conc == 0.25) |
                    (drug == 'Kan' & rel_conc == 0.25) |
                    (drug == 'Neo' & rel_conc == 0.25) |
                    (drug == 'Paro' & rel_conc == 0.125) |
                    (drug == 'Ribo' & rel_conc == 0.125)) %>%
  select(drug, position, wt, mut, score)

# Test correlation between drugs:
dm_data_drug <- pivot_wider(dm_data, names_from = drug, values_from = score)
p_drug_cors <- ggpairs(dm_data_drug, columns = c('Ami', 'G418', 'Kan', 'Neo', 'Paro', 'Ribo'),
                       title = "Correlation between drug scores at best concentrations from Melnikov et al. 2014 (APH(3')-II)")
ggsave('figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_correlation.pdf', p_drug_cors, units = 'cm', width = 25, height = 25)

# All but Ami seem to correlate well -> filter Ami and average rest?
# Could also take worst or something to account
dm_data <- filter(dm_data, !drug == 'Ami') %>%
  group_by(position, wt, mut) %>%
  summarise(score = mean(score, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(raw_score = score,
         class = get_variant_class(wt, mut))

########
