"""
Rules for analysing combined data and calculating AA subtypes

Expects global variables:
"""

#### Analyse whole dataset ####
rule landscape_dimensionality_reduction:
    """
    Analyse the mutational landscape using PCA, tSNE and UMAP dimensionality
    reduction, relating them to various biophysical properties
    """
    input:
         'data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape/pc1_vs_mean_score.pdf',
        'figures/1_landscape/pc1_pc2_mean_score.pdf',
        'figures/1_landscape/pc2_pc4_surface_accessibility.pdf',
        'figures/1_landscape/pc2_vs_surface_accessibility.pdf',
        'figures/1_landscape/pc2_vs_hydrophobicity.pdf',
        'figures/1_landscape/foldx_pc_cor.pdf',
        'figures/1_landscape/tsne_study.pdf',
        'figures/1_landscape/tsne_aa.pdf',
        'figures/1_landscape/tsne_hydrophobicity.pdf',
        'figures/1_landscape/tsne_surface_accessibility.pdf',
        'figures/1_landscape/umap_study.pdf',
        'figures/1_landscape/umap_aa.pdf',
        'figures/1_landscape/umap_hydrophobicity.pdf',
        'figures/1_landscape/umap_mean_er.pdf',
        'figures/1_landscape/umap_surface_accessibility.pdf'

    log:
        'logs/landscape_dimensionality_reduction.log'

    shell:
        'Rscript bin/analysis/1_landscape/dimensionality_reduction.R &> {log}'

#### Clustering ####
cluster_plots = ['ramachanran_angles.pdf', 'cluster_sizes.pdf', 'mean_profiles.pdf',
                 'profile_correlation.pdf', 'foldx_profiles.pdf', 'chem_env_profiles.pdf',
                 'clustering.pdf', 'umap.pdf']
rule make_subtypes:
    """
    Generalised clustering script, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/{name}.yaml'

    output:
        'data/subtypes/{name}.tsv',
        [f'figures/2_subtypes/{{name}}/{x}' for x in cluster_plots]

    log:
        'logs/make_subtypes/{name}.log'

    shell:
        'Rscript bin/analysis/2_subtypes/make_subtypes.R --data data/subtypes --figures figures/2_subtypes {input.yaml} &> {log}'

rule characterise_subtypes:
    """
    Make characterisation plots for each amino acid from an input subtype clustering
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/{name}.tsv"

    output:
        [f"figures/2_subtypes/{{name}}/aa_profiles/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_subtypes/{{name}}/aa_profiles_relative/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_subtypes/{{name}}/ss_probabilities/{x}.pdf" for x in AA_ALPHABET],
        "figures/2_subtypes/{name}/er_vs_surface_accessibility.pdf",
        "figures/2_subtypes/{name}/er_vs_size.pdf"

    log:
        'logs/characterise_subtypes/{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes {input.subtypes} &> {log}"

rule all_position_subtypes:
    """
    Genreate subtypes from positions of all amino acids (using Hclust with dynamic cutting)
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "data/clustering/hclust_profile_dynamic_all_positions.tsv",
        [f"figures/2_subtypes/hclust_profile_dynamic_all_positions/{x}" for x in cluster_plots],
        "figures/2_subtypes/hclust_profile_dynamic_all_positions/cluster_occupancy.pdf"

    log:
        "logs/all_position_subtypes.log"

    shell:
        "Rscript bin/analysis/2_subtypes/all_positions.R &> {log}"

rule all_position_characterisation:
    """
    Analyse all positions against the various biochemical factors directly
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "figures/3_continuous/mean_er_correlations.pdf"

    log:
        "logs/all_position_characterisation.log"

    shell:
        "Rscript bin/analysis/3_continuous/all_position_characterisation.R &> {log}"