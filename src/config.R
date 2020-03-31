#!/usr/bin/env Rscript 
# Load libraries and configurations for all analyses, to keep consistency within project

# Load all libraries (keeps namespace consistent)
library(Biostrings)
library(mclust)

library(tidyverse)
library(rlang)
library(magrittr)
library(broom)
library(readxl)
library(ggpubr)
library(GGally)
library(ggtext)
library(ggdendro)
library(grImport2)
library(grConvert)
library(uwot)
library(yaml)
library(dynamicTreeCut)
library(dbscan)
library(multipanelfigure)
library(tools)
library(cluster)

# Custom packages, available at github.com/allydunham/XXX
library(tblhelpr)
library(plotlistr)

# Universal project functions
source('src/subtypes_utils.R')

# Gene Uniprot IDS
UNIPROT_IDS <- structure(sapply(dir('data/studies/'), function(x){read_yaml(str_c('data/studies/', x, '/', x, '.yaml')) %>% extract2('uniprot_id')}),
                         names = sapply(dir('data/studies/'), function(x){read_yaml(str_c('data/studies/', x, '/', x, '.yaml')) %>% extract2('gene')}))

### GGPlot theme ###
# clean with centered title by default
theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

### Colour Scheme ###
MUT_CLASS_COLOURS <- c(Missense='cornflowerblue', Nonsense='firebrick2', Synonymous='green2')

# AA colours, roughly by chemistry
AA_COLOURS <- c(A='red', I='salmon', L='firebrick', M='orange', V='tomato',
                F='gold', W='yellow3', Y='khaki3',
                N='cadetblue1', C='cornflowerblue', Q='cyan', S='blue', T='darkslateblue',
                R='green', H='green4', K='seagreen1',
                D='purple', E='pink',
                G='antiquewhite2', P='black', X='grey')

# Colourbrewer scales for various metrics
CLUSTER_COLOURS <- list(type='qual', palette='Set3', direction = -1)
ER_PROFILE_COLOURS <- list(type='div', palette='RdBu', direction = 1)
ER_COR_COLOURS <- list(type='div', palette='RdYlBu', direction = 1)
ER_DIST_COLOURS <- list(type='seq', palette='BuPu', direction = 1)
SIFT_COLOURS <- list(type='seq', palette='PuRd', direction = -1)
FOLDX_COLOURS <- list(type='div', palette='PiYG', direction = 1)
SS_COLOURS <- list(type='div', palette='BrBG', direction = 1)
CHEM_ENV_COLOURS <- list(type='seq', palette='GnBu', direction = 1)
AA_DISTANCE_COLOURS <- list(type='seq', palette='YlOrRd', direction = -1)

ER_COSINE_COLOURS <- list(colours = c('#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#000000', 
                                      '#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7'),
                          limits = c(-1.001, 1.001)) # Uses scale_x_gradientn instead of brewer

### Biological Data ###
# Categories of amino acid
# From Sigma-Aldrich website (https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
AA_REDUCED_CLASSES <- list(Aliphatic=c('A', 'I', 'L', 'M', 'V'),
                           Aromatic=c('F', 'W', 'Y'),
                           Polar=c('N', 'C', 'Q', 'S', 'T'),
                           Basic=c('R', 'H', 'K'),
                           Acidic=c('D', 'E'),
                           Glycine='G',
                           Proline='P')
AA_REDUCED_HASH <- structure(rep(names(AA_REDUCED_CLASSES), times=sapply(AA_REDUCED_CLASSES, length)),
                             names=unlist(AA_REDUCED_CLASSES))

### Output Text Encoding ###
# Where applicable encodings have a plaintext and plotmath string version.
# The plotmath strings are converted to expressions wherenecessary with something like sapply(X_PLOTMATH, function(x){parse(text = x)})

# FoldX Terms
FOLDX_TERMS <- c(total_energy='Total Energy', backbone_hbond='Backbone H-bonds', sidechain_hbond='Sidechain H-bonds', van_der_waals='Van der Waals Forces', 
                 electrostatics='Electrostatics', solvation_polar='Solvation (Polar)', solvation_hydrophobic='Solvation (Hydrophobic)', 
                 van_der_waals_clashes='Van der Waals Clashes', entropy_sidechain='Sidechain Entropy', entropy_mainchain='Mainchain Entropy', 
                 cis_bond='cis-Peptide Bonds', torsional_clash='Torsional Clash', backbone_clash='Backbone Clash', helix_dipole='Helix Dipoles',
                 disulfide='Disulphide Bonds', electrostatic_kon='Electorstatic Kon', partial_covalent_bonds='Partical Covalent Bonds', 
                 energy_ionisation='Ionisation Energy', sloop_entropy='S-Loop Entropy', mloop_entropy='M-Loop Entropy', entropy_complex='Complex Entropy',
                 water_bridge='Water Bridges')
FOLDX_TERMS_PLOTMATH <- str_c("'", FOLDX_TERMS, "'") %>% set_names(names(FOLDX_TERMS))
FOLDX_TERMS_PLOTMATH['cis_bond'] <- "italic(cis)*'-Peptide Bonds'"

# DSSP Classes
DSSP_CLASSES <- c(g="3_10 Helix", h="Alpha Helix", i="Pi Helix", e="Beta Strand", b="Beta Bridge", c="Coil", s="Bend", t="Turn")
DSSP_CLASSES_PLOTMATH <- c(g="3[10]*'-helix'", h="alpha*'-helix'", i="pi*'-helix'", e="beta*'-strand'", b="beta*'-bridge'", c="'Coil'", s="'Bend'", t="'Turn'")

# Clusters
# Clusters are named X1, X2, ... for regular subtypes, XO for outliers and XP for designated 'permissive' subtypes
# This scheme is applied after the main algorithm has been run, with permissive positions designated separatly and simple integer encoding (0 = outlier) used
# during the algorithm itself
CLUSTER_OUTLIER_CHAR <- 'O'
CLUSTER_PERMISSIVE_CHAR <- 'P'

CLUSTER_RE <- str_c('[A-Y](?:[0-9]*|', CLUSTER_OUTLIER_CHAR,'|', CLUSTER_PERMISSIVE_CHAR,')')
CLUSTER_OUTLIER_RE <- str_c('[A-Y]', CLUSTER_OUTLIER_CHAR)
CLUSTER_PERMISSIVE_RE <- str_c('[A-Y]', CLUSTER_PERMISSIVE_CHAR)
