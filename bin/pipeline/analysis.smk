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

#### Clustering ####
cluster_plots = ['ramachanran_angles.pdf', 'cluster_sizes.pdf', 'mean_profiles.pdf',
                 'profile_correlation.pdf', 'foldx_profiles.pdf', 'chem_env_profiles.pdf']
rule kmeans_clustering:
    """
    Perform K-means clustering on each amino acid type to break them down into subtypes, and
    produce plots analysing the clustering
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/kmeans_{mode}_k_{k}_min_{min}.tsv',
        [f'figures/2_clustering/kmeans_{{mode}}_k_{{k}}_min_{{min}}/{x}' for x in cluster_plots]

    log:
        'logs/kmeans_clustering/{mode}_k_{k}_min_{min}.log'

    shell:
        'Rscript bin/analysis/2_clustering/kmeans_clustering.R --ncluster {wildcards.k} --min_size {wildcards.min} --mode {wildcards.mode} &> {log}'

rule hclust_clustering:
    """
    Perform Heirarchical clustering on each amino acid type to break them down into subtypes, and
    produce plots analysing the clustering
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/hclust_{mode}_{cut}_{h}_min_{min}_distance_{distance}.tsv',
        [f'figures/2_clustering/hclust_{{mode}}_{{cut}}_{{h}}_min_{{min}}_distance_{{distance}}/{x}' for x in cluster_plots]

    log:
        'logs/hclust_clustering/{mode}_{cut}_{h}_min_{min}_distance_{distance}.log'

    shell:
        'Rscript bin/analysis/2_clustering/hclust_clustering.R --{wildcards.cut} {wildcards.h} --min_size {wildcards.min} --distance {wildcards.distance} --mode {wildcards.mode} &> {log}'

rule hdbscan_clustering:
    """
    Perform HDBSCAN clustering on each amino acid type to break them down into subtypes, and
    produce plots analysing the clustering
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/hdbscan_{mode}_min_{min}_distance_{distance}.tsv',
        [f'figures/2_clustering/hdbscan_{{mode}}_min_{{min}}_distance_{{distance}}/{x}' for x in cluster_plots]

    log:
        'logs/hdbscan_clustering/{mode}_min_{min}_distance_{distance}.log'

    shell:
        'Rscript bin/analysis/2_clustering/hdbscan_clustering.R --min_size {wildcards.min} --mode {wildcards.mode} --distance {wildcards.distance} &> {log}'

rule dbscan_clustering:
    """
    Perform DBSCAN clustering on each amino acid type to break them down into subtypes, and
    produce plots analysing the clustering
    """
    input:
        'data/combined_mutational_scans.tsv'

    output:
        'data/clusterings/dbscan_{mode}_min_{min}_eps_{eps}_distance_{distance}.tsv',
        [f'figures/2_clustering/dbscan_{{mode}}_min_{{min}}_eps_{{eps}}_distance_{{distance}}/{x}' for x in cluster_plots]

    log:
        'logs/dbscan_clustering/{mode}_min_{min}_eps_{eps}_distance_{distance}.log'

    shell:
        'Rscript bin/analysis/2_clustering/dbscan_clustering.R --minPts {wildcards.min} --eps {wildcards.eps} --mode {wildcards.mode} --distance {wildcards.distance} &> {log}'
