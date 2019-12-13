"""
Rules for the study standardisation pipeline and generrating summary information on the data

Expects global variables:
  - STUDIES: dict of study properties
  - GENES: dict linking genes to lists of studies
  - UNFILTERED_STUDIES: list of unfiltered studies
  - UNFILTERED_GENES: list of unfiltered genes
"""

#### Summarisation ####
rule summarise_study_set:
    """
    Produce tables detailing the studies processed for the project
    """
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        study='meta/study_summary.tsv',
        gene='meta/gene_summary.tsv',
        overall='meta/overall_summary',

    log:
        "logs/summarise_study_set.log"

    shell:
        "python bin/analysis/0_data_properties/summarise_studies.py -s {output.study} -g {output.gene} -u {output.overall} data/studies/* &> {log}"

rule study_summary_plots:
    """
    Make plots summarising the dataset
    """
    input:
        'meta/study_summary.tsv',
        'meta/gene_summary.tsv',
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        'figures/0_data_properties/study_variants_summary.pdf',
        'figures/0_data_properties/gene_variants_summary.pdf',
        'figures/0_data_properties/position_coverage.pdf'

    log:
        "logs/study_summary_plots.log"

    shell:
        'Rscript bin/analysis/0_data_properties/data_summary_plots.R &> {log}'

rule summarise_standardised_data:
    """
    Make plots summarising the standardised and filtered dataset
    """
    input:
        'data/long_combined_mutational_scans.tsv'

    output:
        'figures/0_data_properties/standardised_distributions.pdf',
        'figures/0_data_properties/position_data_summary.pdf'

    log:
        "logs/summarise_standardised_data.log"

    shell:
        'Rscript bin/analysis/0_data_properties/summarise_standardised_data.R &> {log}'


#### Combine Deep Mutational Scans ####
rule standardise_study:
    """
    Import a study and convert it into the standardised table format
    """
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
        "Rscript {input} &> {log}"

rule standardise_all_studies:
    """
    Run standardisation procedure for all studies
    """
    input:
        expand('data/studies/{study}/{study}.tsv', study=STUDIES.keys())

rule combine_dms_data:
    """
    Gather all standardised data and additional data on the genes and structures and
    combine it into a single tsv, including calculating principle components and tSNE
    projection for the wide version of the data.
    """
    input:
        expand('data/studies/{study}/{study}.{ext}', study=UNFILTERED_STUDIES, ext=('tsv', 'yaml')),
        expand('data/fasta/{gene}.fa', gene=UNFILTERED_GENES),
        expand('data/sift/{gene}.SIFTprediction', gene=UNFILTERED_GENES),
        expand('data/foldx/{gene}/average_{gene}.fxout', gene=UNFILTERED_GENES),
        expand('data/backbone_angles/{gene}.tsv', gene=UNFILTERED_GENES),
        expand('data/surface_accessibility/{gene}.rsa', gene=UNFILTERED_GENES),
        expand('data/chemical_environment/{gene}_within_10.0.tsv', gene=UNFILTERED_GENES),
        expand('data/porter5/{gene}.ss8', gene=UNFILTERED_GENES),
        'meta/residue_hydrophobicity.tsv'

    output:
        'data/long_combined_mutational_scans.tsv',
        'data/combined_mutational_scans.tsv'

    log:
        "logs/combine_dms_data.log"

    shell:
        f"Rscript bin/data_processing/combine_standardised_data.R {' '.join([f'data/studies/{s}' for s in UNFILTERED_STUDIES])} &> {{log}}"
