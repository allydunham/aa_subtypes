"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""
# General todo list for pipeline
# TODO Better logging
# TODO Better all rules
# TODO rename 1_dimensionality_reduction?
# TODO add automated download of some of the input data?
# TODO add docstrings to rules

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
    make_gene_fasta, all_foldx_predictions, foldx_variants, foldx_split

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
include: 'bin/pipeline/sift.smk'
include: 'bin/pipeline/foldx.smk'
include: 'bin/pipeline/analysis.smk'

#### Global rules ####
# TODO - group up plots into lists?
rule all:
    input:
        'data/combined_mutational_scans.tsv' # Covers all standardisation, SIFT and FoldX
        'meta/study_summary.tsv',
        'meta/gene_summary.tsv',
        'meta/overall_summary',
        VALIDATION_PLOTS,
        'figures/0_data_properties/study_variants_summary.pdf',
        'figures/0_data_properties/gene_variants_summary.pdf',
        'figures/0_data_properties/position_coverage.pdf',
        'data/clusterings/kmeans_3.tsv'

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

#### Summarise dataset ####
rule summarise_dataset:
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        study='meta/study_summary.tsv',
        gene='meta/gene_summary.tsv',
        overall='meta/overall_summary',

    shell:
        "python bin/utils/summarise_studies.py -s {output.study} -g {output.gene} -u {output.overall} data/studies/*"

rule data_summary_plots:
    input:
        'meta/study_summary.tsv',
        'meta/gene_summary.tsv',
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        'figures/0_data_properties/study_variants_summary.pdf',
        'figures/0_data_properties/gene_variants_summary.pdf',
        'figures/0_data_properties/position_coverage.pdf'

    shell:
        'Rscript bin/analysis/0_data_properties/data_summary_plots.R'

#### Combine Deep Mutational Scans ####
rule standardise_study:
    input:
        "data/studies/{study}/standardise_{study}.R",
        lambda wildcards: [f'data/studies/{wildcards.study}/raw/{x}' for x in
                           STUDIES[wildcards.study]['input_files']]

    output:
        "data/studies/{study}/{study}.tsv",
        "figures/0_data_properties/per_study/{study}/original_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/transformed_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/normalised_distribution.pdf"

    log:
        "logs/standardise_study/{study}.log"

    shell:
        "Rscript {input} 2> {log}"

rule combine_dms_data:
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('tsv', 'yaml')),
        expand('data/sift/{gene}.{ext}', gene=GENES.keys(), ext=('fa', 'SIFTPrediction')),
        expand('data/foldx/{gene}/average_{gene}.fxout', gene=GENES.keys()),
        expand('data/backbone_angles/{gene}.tsv', gene=GENES.keys())

    output:
        'data/combined_mutational_scans.tsv'

    shell:
        f"Rscript bin/data_processing/combine_standardised_data.R {' '.join([f'data/studies/{s}' for s,v in STUDIES.items() if not v['qc']['filter']])}"

#### Misc Property calculations ####
rule calculate_backbone_angles:
    input:
        pdb="data/foldx/{gene}/{gene}_Repair.pdb",
        yaml="meta/structures.yaml"

    output:
        "data/backbone_angles/{gene}.tsv"

    log:
        "logs/calculate_backbone_angles/{gene}.log"

    shell:
        "python bin/data_processing/get_backbone_angles.py --yaml {input.yaml} {input.pdb} > {output} 2> {log}"

rule filter_pdb:
    input:
        'data/foldx/{gene}/{gene}_Repair.pdb'

    output:
        'data/surface_accessibility/{gene}.pdb'

    log:
        'logs/filter_pdb/{gene}.log'

    shell:
        'python bin/data_processing/filter_pdb.py --yaml meta/structures.yaml {input} > {output} 2> {log}'

rule naccess:
    input:
        rules.filter_pdb.output

    output:
        asa='data/surface_accessibility/{gene}.asa',
        rsa='data/surface_accessibility/{gene}.rsa'

    log:
        'logs/naccess/{gene}.log'

    shell:
        """
        naccess {input}
        mv {wildcards.gene}.log {log}
        mv {wildcards.gene}.asa {output.asa}
        mv {wildcards.gene}.rsa {output.rsa}
        """