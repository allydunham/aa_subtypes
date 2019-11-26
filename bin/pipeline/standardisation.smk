"""
Rules for the study standardisation pipeline and generrating summary information on the data
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

    shell:
        "python bin/utils/summarise_studies.py -s {output.study} -g {output.gene} -u {output.overall} data/studies/*"

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

    shell:
        'Rscript bin/analysis/0_data_properties/data_summary_plots.R'

rule summarise_standardised_data:
    """
    Make plots summarising the standardised and filtered dataset
    """
    input:
        'data/long_combined_mutational_scans.tsv'

    output:
        'figures/0_data_properties/standardised_distributions.pdf',
        'figures/0_data_properties/position_data_summary.pdf'

    shell:
        'Rscript bin/analysis/0_data_properties/summarise_standardised_data.R'


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
        "Rscript {input} 2> {log}"

rule combine_dms_data:
    """
    Gather all standardised data and additional data on the genes and structures and
    combine it into a single tsv, including calculating principle components and tSNE
    projection for the wide version of the data.
    """
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('tsv', 'yaml')),
        expand('data/sift/{gene}.{ext}', gene=GENES.keys(), ext=('fa', 'SIFTprediction')),
        expand('data/foldx/{gene}/average_{gene}.fxout', gene=GENES.keys()),
        expand('data/backbone_angles/{gene}.tsv', gene=GENES.keys()),
        expand('data/surface_accessibility/{gene}.rsa', gene=GENES.keys()),
        expand('data/chemical_environment/{gene}_within_10.0.tsv', gene=GENES.keys()),
        'meta/residue_hydrophobicity.tsv'

    output:
        'data/long_combined_mutational_scans.tsv',
        'data/combined_mutational_scans.tsv'

    shell:
        f"Rscript bin/data_processing/combine_standardised_data.R {' '.join([f'data/studies/{s}' for s,v in STUDIES.items() if not v['qc']['filter']])}"
