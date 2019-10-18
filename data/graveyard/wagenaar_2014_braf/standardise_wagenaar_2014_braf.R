#!/usr/bin/env Rscript
# Standardise data from Wagenaar et al. 2014 (BRAF)
source('src/config.R')
source('src/study_standardising.R')

# Import and process data
meta <- read_yaml('data/studies/wagenaar_2014_braf/wagenaar_2014_braf.yaml')
dm_data <- read_xls('data/studies/wagenaar_2014_braf/raw/wagenaar_2014_braf.xls', skip = 3) %>%
  rename(position = Position,
         mut = acid,
         median_enrichment = Median,
         rep1_codon1 = `Replicate 1`,
         rep1_codon2 = ...5,
         rep1_codon3 = ...6,
         rep1_codon4 = ...7,
         rep1_codon5 = ...8,
         rep1_codon6 = ...9,
         rep2_codon1 = `Replicate 2`,
         rep2_codon2 = ...11,
         rep2_codon3 = ...12,
         rep2_codon4 = ...13,
         rep2_codon5 = ...14,
         rep2_codon6 = ...15,
         ic50_vs_brafV600E = BRAFV600E,
         individually_tested = `mutant?`,
         possible_by_single_sub = `substitution?`) %>%
  filter(!is.na(rep1_codon1) & !rep1_codon1 == 'Replicate 1') %>%
  mutate_at(vars(-mut, -individually_tested, -possible_by_single_sub, -ic50_vs_brafV600E), as.numeric)%>%
  mutate(wt = str_split(meta$seq, '')[[1]][position],
         score = log2(median_enrichment),
         raw_score = median_enrichment,
         class = get_variant_class(wt, mut)) %>%
  select(position, wt, mut, score, raw_score, class)

# Save output
standardise_study(dm_data, meta$study, meta$transform)