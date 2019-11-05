# pylint: disable-all
"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""
import os
from Bio import SeqIO
from ruamel.yaml import YAML

import subtypes_utils as sutil

yaml = YAML(typ='safe')
UNIREF90_DB_PATH = '/hps/research1/beltrao/ally/databases/uniref90/uniref90_2019_1.fasta'
FASTA_LINE_LENGTH = 80

configfile: 'snakemake.yaml'
localrules: all, all_standardisation, all_sift, make_sift_fastas

rule all:
    input:
        expand('data/studies/{study}/{study}.tsv', study=config['studies']),
        expand('data/sift/{gene}.SIFTprediction', gene=config['genes'].keys())

rule all_standardisation:
    input:
        expand('data/studies/{study}/{study}.tsv', study=config['studies'])

rule all_sift:
    input:
        expand('data/sift/{gene}.SIFTprediction', gene=config['genes'].keys())

#### Validate Data ####
# Validate Melnikov et al. 2014 (APH(3')-II)
rule validate_melnikov:
    input:
        config['input_files']['melnikov_2014_aph3ii']

    output:
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/initial_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/filtered_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/rel_conc_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_correlation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_melnikov_2014_aph3ii.R"

# Validate Kitzman et al. 2015 (GAL4)
rule validate_kitzman:
    input:
        config['input_files']['kitzman_2015_gal4']

    output:
        "figures/0_data_properties/per_study/kitzman_2015_gal4/validate_selection_combination.pdf"

    script:
        "bin/analysis/0_data_properties/validate_kitzman_2015_gal4.R"

# Validate Giacomelli et al. 2018 (TP53)
rule validate_giacomelli:
    input:
        config['input_files']['giacomelli_2018_tp53']

    output:
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/initial_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/conditions.pdf"

    script:
        "bin/analysis/0_data_properties/validate_giacomelli_2018_tp53.R"

# Validate Heredia et al. 2018
rule validate_heredia:
    input:
       config['input_files']['heredia_2018_ccr5'] + config['input_files']['heredia_2018_cxcr4']

    output:
        "figures/0_data_properties/per_study/heredia_2018_ccr5/replicate_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_ccr5/experiment_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_cxcr4/replicate_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_cxcr4/experiment_correlation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_heredia_2018.R"

# Validate Sarkisyan et al. 2016 (GFP)
rule validate_sarkisyan:
    input:
        config['input_files']['sarkisyan_2016_gfp']

    output:
        "figures/0_data_properties/per_study/sarkisyan_2016_gfp/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_sarkisyan_2016_gfp.R"

# Validate Dorrity et al. 2018 (STE12)
rule validate_dorrity:
    input:
        config['input_files']['dorrity_2018_ste12']

    output:
        "figures/0_data_properties/per_study/dorrity_2018_ste12/rep_correlation.pdf",
        "figures/0_data_properties/per_study/dorrity_2018_ste12/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_dorrity_2018_ste12.R"

# Validate Araya et al. 2012 (YAP1)
rule validate_araya:
    input:
        config['input_files']['araya_2012_yap1']

    output:
        "figures/0_data_properties/per_study/araya_2012_yap1/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_araya_2012_yap1.R"

# Validate Starita et al. 2013 (UBE4B)
rule validate_starita:
    input:
        config['input_files']['starita_2013_ube4b']

    output:
        "figures/0_data_properties/per_study/starita_2013_ube4b/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_starita_2013_ube4b.R"

# Check correlation with SIFT
rule sift_correlation:
    input:
        expand('data/studies/{study}/{study}.tsv', study=config['studies']),
        expand('data/studies/{study}/{study}.yaml', study=config['studies']),
        expand('data/sift/{gene}.SIFTprediction', gene=config['genes'].keys()),
        expand('data/sift/{gene}.fa', gene=config['genes'].keys())

    output:
        'figures/0_data_properties/sift_score_correlation.pdf'

    script:
        'bin/analysis/0_data_properties/sift_correlation.R'

#### Standardise Data ####
# Process the raw data from each study
rule standardise_study:
    input:
        "data/studies/{study}/standardise_{study}.R",
        lambda wildcards: [f'data/studies/{wildcards.study}/raw/{x}' for x in
                           config['input_files'][wildcards.study]]

    output:
        "data/studies/{study}/{study}.tsv",
        "figures/0_data_properties/per_study/{study}/original_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/transformed_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/normalised_distribution.pdf"

    log:
        "logs/standardise_study/{study}.log"

    shell:
        "Rscript {input} 2> {log}"

#### Make Tool Predictions ####
rule make_gene_fasta:
    input:
        lambda wildcards: [f'data/studies/{i}/{i}.yaml' for i in config['genes'][wildcards.gene]]

    output:
        "data/sift/{gene}.fa"

    run:
        seq = None
        for study_yaml in input:
            with open(study_yaml, 'r') as yaml_file:
                conf = yaml.load(yaml_file)

            if seq is None:
                seq = conf['seq']
            elif not seq == conf['seq']:
                raise ValueError(f"Two studies on {wildcards.gene} have different sequences")

        with open(output[0], 'w') as fasta_file:
            print(f">{wildcards.gene}", file=fasta_file)
            for i in range(0, len(seq), FASTA_LINE_LENGTH):
                print(seq[i:(i + FASTA_LINE_LENGTH)], file=fasta_file)

rule sift4g:
    input:
        fa = "data/sift/{gene}.fa",
        db = UNIREF90_DB_PATH

    output:
        "data/sift/{gene}.SIFTprediction"

    log:
        'logs/sift4g/{gene}.log'

    resources:
        mem_mb = 8000

    shell:
        "sift4g -q {input.fa} -d {input.db} --out data/sift 2> {log}"

# Make all FoldX for study genes
# TODO

# Plots