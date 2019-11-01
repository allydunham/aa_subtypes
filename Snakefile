# pylint: disable-all
"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""

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
        expand('data/sift/{gene}.SIFTprediction', gene=config['genes'])

rule all_standardisation:
    input:
        expand('data/studies/{study}/{study}.tsv', study=config['studies'])

rule all_sift:
    input:
        expand('data/sift/{gene}.SIFTprediction', gene=config['genes'])

#### Validate Data ####
# Test multiple mutation averaging
rule validate_multi_muts:
    input:
        "data/studies/starita_2013_ube4b/raw/starita_2013_ube4b_ubox.xlsx",
        "data/studies/araya_2012_yap1/raw/araya_2012_hYAP65_ww.tsv"

    output:
        "figures/0_data_properties/averaging_multi_mutants.pdf"

    script:
        "bin/analysis/0_data_properties/validate_multi_muts.R"

# Validate Melnikov et al. 2014 (APH(3')-II)
rule validate_melnikov:
    # Requires Melnikov .aacount files to be in data/studies/melnikov_2014_aph3ii/raw
    output:
        "figures/0_data_properties/melnikov_2014_aph3ii/initial_library_correlation.pdf",
        "figures/0_data_properties/melnikov_2014_aph3ii/filtered_library_correlation.pdf",
        "figures/0_data_properties/melnikov_2014_aph3ii/rel_conc_correlation.pdf",
        "figures/0_data_properties/melnikov_2014_aph3ii/drug_correlation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_melnikov_2014_aph3ii.R"

# Validate Kitzman et al. 2015 (GAL4)
rule validate_kitzman:
    input:
        "data/studies/kitzman_2015_gal4/raw/kitzman_2015_gal4_enrichment.xlsx"

    output:
        "figures/0_data_properties/kitzman_2015_gal4/validate_selection_combination.pdf"

    script:
        "bin/analysis/0_data_properties/validate_kitzman_2015_gal4.R"

# Validate Giacomelli et al. 2018 (TP53)
rule validate_giacomelli:
    input:
        "data/studies/giacomelli_2018_tp53/raw/41588_2018_204_MOESM5_ESM.xlsx"

    output:
        "figures/0_data_properties/giacomelli_2018_tp53/initial_experiment_cor.pdf",
        "figures/0_data_properties/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf",
        "figures/0_data_properties/giacomelli_2018_tp53/conditions.pdf"

    script:
        "bin/analysis/0_data_properties/validate_giacomelli_2018_tp53.R"

# Validate Heredia et al. 2018
rule validate_heredia:
    input:
        "data/studies/heredia_2018_ccr5/raw/GSE100368_enrichment_ratios_CCR5.xlsx"

    output:
        "figures/0_data_properties/heredia_2018_ccr5/replicate_correlation.pdf",
        "figures/0_data_properties/heredia_2018_ccr5/experiment_correlation.pdf",
        "figures/0_data_properties/heredia_2018_cxcr4/replicate_correlation.pdf",
        "figures/0_data_properties/heredia_2018_cxcr4/experiment_correlation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_heredia_2018.R"

# Validate Sarkisyan et al. 2016 (GFP)
rule validate_sarkisyan:
    input:
        "data/studies/sarkisyan_2016_gfp/raw/sarkisyan_2016_gfp_AAs.tsv"

    output:
        "figures/0_data_properties/sarkisyan_2016_gfp/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_sarkisyan_2016_gfp.R"

# Validate Dorrity et al. 2018 (STE12)
rule validate_dorrity:
    input:
        "data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd01.xlsx",
        "data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd02.xlsx"

    output:
        "figures/0_data_properties/dorrity_2018_ste12/rep_correlation.pdf",
        "figures/0_data_properties/dorrity_2018_ste12/multi_mut_validation.pdf"

    script:
        "bin/analysis/0_data_properties/validate_dorrity_2018_ste12.R"

#### Standardise Data ####
# Process the raw data from each study
rule standardise_study:
    input:
        "data/studies/{study}/standardise_{study}.R"

    output:
        "data/studies/{study}/{study}.tsv",
        "figures/0_data_properties/{study}/original_distribution.pdf",
        "figures/0_data_properties/{study}/transformed_distribution.pdf",
        "figures/0_data_properties/{study}/normalised_distribution.pdf"

    log:
        "logs/standardise_study/{study}.log"

    shell:
        "Rscript {input} 2> {log}"

#### Make Tool Predictions ####
# Make all SIFT predictions for study genes
# TODO currently overwrites all and re-runs sift every time
rule make_sift_fastas:
    input:
        expand('data/studies/{study}/{study}.yaml', study=config['studies'])

    output:
        expand("data/sift/{gene}.fa", gene=config['genes'])

    run:
        genes = []
        for study_yaml in input:
            with open(study_yaml, 'r') as yaml_file:
                conf = yaml.load(yaml_file)

            if conf['gene'] in genes:
                continue

            genes.append(conf['gene'])
            gene_filesafe = sutil.gene_to_filename(conf['gene'])
            with open(f"data/sift/{gene_filesafe}.fa", 'w') as fasta_file:
                print(f">{gene_filesafe}", file=fasta_file)
                for i in range(0, len(conf['seq']), FASTA_LINE_LENGTH):
                    print(conf['seq'][i:(i + FASTA_LINE_LENGTH)], file=fasta_file)

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