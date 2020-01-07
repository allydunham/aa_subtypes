"""
Rules for analysing combined data and calculating AA subtypes

Expects global variables:
"""

#### Analyse whole dataset ####
rule principle_component_analysis:
    """
    Analyse principle components of the mutational landscape, relating them
    to various biophysical properties
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape_properties/pc1_vs_mean_score.pdf',
        'figures/1_landscape_properties/pc1_pc2_mean_score.pdf',
        'figures/1_landscape_properties/pc2_pc4_surface_accessibility.pdf',
        'figures/1_landscape_properties/pc2_vs_surface_accessibility.pdf',
        'figures/1_landscape_properties/pc2_vs_hydrophobicity.pdf',
        'figures/1_landscape_properties/foldx_pc_cor.pdf'

    log:
        'logs/principle_component_analysis.log'

    shell:
        'Rscript bin/analysis/1_landscape_properties/principle_component_analysis.R &> {log}'

rule tsne_analysis:
    """
    Analyse tSNE dimensionality reduction of the mutational landscape, relating it to biophysical
    properties
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape_properties/tsne_study.pdf',
        'figures/1_landscape_properties/tsne_aa.pdf',
        'figures/1_landscape_properties/tsne_hydrophobicity.pdf',
        'figures/1_landscape_properties/tsne_surface_accessibility.pdf'

    log:
        'logs/tsne_analysis.log'

    shell:
        'Rscript bin/analysis/1_landscape_properties/tsne.R &> {log}'

rule umap_analysis:
    """
    Analyse UMAP dimensionality reduction of the mutational landscape, relating it to biophysical
    properties
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape_properties/umap_study.pdf',
        'figures/1_landscape_properties/umap_aa.pdf',
        'figures/1_landscape_properties/umap_hydrophobicity.pdf',
        'figures/1_landscape_properties/umap_mean_er.pdf',
        'figures/1_landscape_properties/umap_surface_accessibility.pdf'

    log:
        'logs/umap_analysis.log'

    shell:
        'Rscript bin/analysis/1_landscape_properties/umap.R &> {log}'


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
        yaml='meta/clustering/{name}.yaml'

    output:
        'data/clustering/{name}.tsv',
        [f'figures/2_clustering/{{name}}/{x}' for x in cluster_plots]

    log:
        'logs/make_subtypes/{name}.log'

    shell:
        'Rscript bin/analysis/2_clustering/make_subtypes.R --data data/clustering --figures figures/2_clustering {input.yaml} &> {log}'

rule characterise_subtypes:
    """
    Make characterisation plots for each amino acid from an input subtype clustering
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/clustering/{name}.tsv"

    output:
        [f"figures/2_clustering/{{name}}/aa_profiles/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_clustering/{{name}}/aa_profiles_relative/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_clustering/{{name}}/ss_probabilities/{x}.pdf" for x in AA_ALPHABET],
        "figures/2_clustering/{name}/er_vs_surface_accessibility.pdf",
        "figures/2_clustering/{name}/er_vs_size.pdf"

    log:
        'logs/characterise_subtypes/{name}.log'

    shell:
        "Rscript bin/analysis/2_clustering/characterise_subtypes.R --dms {input.dms} --figures figures/2_clustering {input.subtypes} &> {log}"

rule all_position_subtypes:
    """
    Genreate subtypes from positions of all amino acids (using Hclust with dynamic cutting)
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "data/clustering/hclust_profile_dynamic_all_positions.tsv",
        [f"figures/2_clustering/hclust_profile_dynamic_all_positions/{x}" for x in cluster_plots],
        "figures/2_clustering/hclust_profile_dynamic_all_positions/cluster_occupancy.pdf"

    log:
        "logs/all_position_subtypes.log"

    shell:
        "Rscript bin/analysis/2_clustering/all_positions.R &> {log}"