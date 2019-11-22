"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""
# General todo list for pipeline
# TODO Better logging
# TODO Better all rules
# TODO add automated download of some of the input data?
# TODO add docstrings to rules
# TODO clean up master pdb file locations?
# TODO general clean up overhaul/check all in order
# TODO rule to setup logging directories?
# TODO Change to only unfiltered gene input in combine_dms_data
# TODO Split out structures.yaml section validation from get_* scripts to utils

import os
import math
from collections import defaultdict

from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from ruamel.yaml import YAML
yaml = YAML(typ='safe')

import subtypes_utils as sutil

configfile: 'snakemake.yaml'
localrules:
    all, quick_clean, full_clean, standardise_all_studies, all_sift_predictions,
    make_gene_fasta, all_foldx_predictions, foldx_variants, foldx_split,
    filter_pdb, naccess

# Hash of study IDs to their config
STUDIES = {}
for study in os.listdir('data/studies'):
    with open(f'data/studies/{study}/{study}.yaml') as yaml_file:
        STUDIES[study] = yaml.load(yaml_file)

# Hash linking genes to studies
GENES = defaultdict(list)
for study, conf in STUDIES.items():
    GENES[sutil.gene_to_filename(conf['gene'])].append(study)

with open('meta/structures.yaml', 'r') as yaml_file:
    STRUCTURES = yaml.load(yaml_file)

AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'

#### Include subroutines ####
include: 'bin/pipeline/data_validation.smk'
include: 'bin/pipeline/standardisation.smk'
include: 'bin/pipeline/sift.smk'
include: 'bin/pipeline/foldx.smk'
include: 'bin/pipeline/structure_statistics.smk'
include: 'bin/pipeline/analysis.smk'

#### Global rules ####
# TODO - group up plots into lists?
rule all:
    input:
        'data/combined_mutational_scans.tsv', # Covers all standardisation, SIFT and FoldX
        'data/long_combined_mutational_scans.tsv',
        rules.summarise_study_set.output,
        VALIDATION_PLOTS,
        rules.study_summary_plots.output,
        rules.summarise_standardised_data.output,
        rules.principle_component_analysis.output,
        rules.tsne_analysis.output,
        'data/clusterings/kmeans_profile_k_4_min_3.tsv',
        'data/clusterings/kmeans_pca_k_4_min_3.tsv',
        'data/clusterings/hclust_profile_height_17_min_3_distance_manhattan.tsv',
        'data/clusterings/hclust_pca_height_6_min_3_distance_manhattan.tsv',
        'data/clusterings/hclust_profile_number_5_min_3_distance_manhattan.tsv',
        'data/clusterings/hclust_pca_number_5_min_3_distance_manhattan.tsv',
        'data/clusterings/hdbscan_profile_min_6_distance_manhattan.tsv',
        'data/clusterings/hdbscan_pca_min_5_distance_manhattan.tsv',
        'data/clusterings/dbscan_profile_min_5_eps_4_distance_manhattan.tsv',
        'data/clusterings/dbscan_pca_min_5_eps_4_distance_manhattan.tsv'

# Only remove rapidly generated results
def quick_clean_files():
    output_files = [f'data/studies/{s}/{s}.tsv' for s in STUDIES.keys()]
    output_files.append('-r figures/0_data_properties/*')
    output_files.append('-r figures/1_landscape_properties/*')
    output_files.append('-r figures/2_clustering/*')
    output_files.append('meta/study_summary.tsv')
    output_files.append('meta/gene_summary.tsv')
    output_files.append('meta/overall_summary')
    output_files.append('data/combined_mutational_scans.tsv')
    output_files.append('data/clustering/*')
    output_files.append('logs/*/*')
    return output_files

rule quick_clean:
    run:
        for i in quick_clean_files():
            shell(f'rm {i} && echo "rm {i}" || true')

# Same as clean, but also remove slower to genreate Sift & FoldX results
rule full_clean:
    run:
        output_files = quick_clean_files()

        # SIFT results
        output_files.extend([f'data/sift/{g}.fa' for g in GENES.keys()])
        output_files.extend([f'data/sift/{g}.SIFTPrediction' for g in GENES.keys()])

        # FoldX results
        output_files.extend([f"data/foldx/{g}/{g}_Repair.pdb" for g in GENES.keys()])
        output_files.extend([f"data/foldx/{g}/individual_list" for g in GENES.keys()])
        output_files.extend([f"data/foldx/{g}/*.fxout" for g in GENES.keys()])
        output_files.extend([f"-r data/foldx/{g}/processing" for g in GENES.keys()])

        for i in output_files:
            shell(f'rm {i} && echo "rm {i}" || true')

rule standardise_all_studies:
    input:
        expand('data/studies/{study}/{study}.tsv', study=STUDIES.keys())

rule all_sift_predictions:
    input:
        expand('data/sift/{gene}.SIFTprediction', gene=GENES.keys())

rule all_foldx_predictions:
    input:
        expand("data/foldx/{gene}/average_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/dif_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/raw_{gene}.fxout", gene=GENES.keys())

rule validate_data:
    input:
        VALIDATION_PLOTS
