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
cluster_plots = ['ramachanran_angles.pdf', 'cluster_sizes.pdf', 'mean_profiles.pdf',
                 'profile_correlation.pdf', 'foldx_profiles.pdf', 'chem_env_profiles.pdf']
rule kmeans_clustering:
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/kmeans_{mode}_{n}.tsv',
        [f'figures/2_clustering/kmeans_{{mode}}_{{n}}/{x}' for x in cluster_plots]

    log:
        'logs/kmeans_clustering/{mode}_{n}.log'

    shell:
        'Rscript bin/analysis/2_clustering/kmeans_clustering.R --ncluster {wildcards.n} --mode {wildcards.mode} &> {log}'

rule hclust_clustering:
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/hclust_{mode}_h_{h}_min_{min}.tsv',
        [f'figures/2_clustering/hclust_{{mode}}_h_{{h}}_min_{{min}}/{x}' for x in cluster_plots]

    log:
        'logs/hclust_clustering/{mode}_h_{h}_min_{min}.log'

    shell:
        'Rscript bin/analysis/2_clustering/hclust_clustering.R --height {wildcards.h} --min_size {wildcards.min} --mode {wildcards.mode} &> {log}'

rule hdbscan_clustering:
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/hdbscan_{mode}_{n}.tsv',
        [f'figures/2_clustering/hdbscan_{{mode}}_{{n}}/{x}' for x in cluster_plots]

    log:
        'logs/hdbscan_clustering/{mode}_{n}.log'

    shell:
        'Rscript bin/analysis/2_clustering/hdbscan_clustering.R --min_size {wildcards.n} --mode {wildcards.mode} &> {log}'
