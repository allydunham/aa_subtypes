"""
Rules for analysing combined data and calculating AA subtypes

Expects global variables:
"""

#### Analyse whole dataset ####
rule principle_component_analysis:
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape_properties/pc1_vs_mean_score.pdf',
        'figures/1_landscape_properties/pc1_pc2_mean_score.pdf',
        'figures/1_landscape_properties/pc2_pc4_surface_accessibility.pdf',
        'figures/1_landscape_properties/pc2_vs_surface_accessibility.pdf',
        'figures/1_landscape_properties/pc2_vs_hydrophobicity.pdf',
        'figures/1_landscape_properties/foldx_pc_cor.pdf'

    shell:
        'Rscript bin/analysis/1_landscape_properties/principle_component_analysis.R'

#### Clustering ####
rule kmeans_clustering:
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/kmeans_{n}.tsv',
        'figures/2_clustering/kmeans_{n}/ramachanran_angles.pdf',
        'figures/2_clustering/kmeans_{n}/cluster_sizes.pdf',
        'figures/2_clustering/kmeans_{n}/mean_profiles.pdf',
        'figures/2_clustering/kmeans_{n}/profile_correlation.pdf'

    log:
        'logs/kmeans_clustering/{n}.log'

    shell:
        'Rscript bin/analysis/2_clustering/kmeans_clustering.R --ncluster {wildcards.n} &> {log}'