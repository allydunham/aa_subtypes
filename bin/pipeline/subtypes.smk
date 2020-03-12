"""
Rules for calculating and analysing AA subtypes

Expects global variables:
- STANDARD_CLUSTERINGS: list of names for clusterings generated
"""

# Standard subtypes plots
diagnostic_plots = ['umap.pdf', 'tsne.pdf', 'silhouette_global.pdf',
                    'silhouette_clustering_vars.pdf', 'silhouette_per_aa.pdf',
                    'silhouette_per_aa_cosine.pdf']

characterisation_plots = ['ramachandran_angles.pdf', 'sizes.pdf', 'er_profiles.pdf',
                          'er_correlation.pdf', 'er_correlation_selective.pdf',
                          'foldx.pdf', 'chem_env.pdf', 'ss_probability.pdf',
                          'er_cosine.pdf', 'er_cosine_selective.pdf',
                          'er_vs_surface_accessibility.pdf', 'er_vs_size.pdf', 'overall_dend.pdf', 'overall_dend_selective.pdf',
                          'multi_position_subtype_consistency.pdf']
characterisation_plots.extend([f"aa_profiles/{x}.pdf" for x in AA_ALPHABET])
characterisation_plots.extend([f"aa_profiles_relative/{x}.pdf" for x in AA_ALPHABET])
characterisation_plots.extend([f"ss_probabilities/{x}.pdf" for x in AA_ALPHABET])
characterisation_plots.extend([f"aa_dends/{x}.pdf" for x in AA_ALPHABET])
characterisation_plots.extend([f"profile_variance/{x}.pdf" for x in AA_ALPHABET])

#### General Rules ####
rule all_position_subtypes:
    """
    Genreate subtypes from positions of all amino acids (using Hclust with dynamic cutting)
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "data/subtypes/all_positions.tsv",
        "data/subtypes/all_positions.rds",
        "figures/2_subtypes/all_positions/umap.pdf",
        "figures/2_subtypes/all_positions/tsne.pdf",
        "figures/2_subtypes/all_positions/silhouette.pdf",
        "figures/2_subtypes/all_positions/cluster_occupancy.pdf"

    log:
        "logs/all_position_subtypes.log"

    shell:
        "Rscript bin/analysis/2_subtypes/all_positions.R &> {log}"

rule evaluate_kmeans_k:
    """
    Assess different values of k to use for k-means clustering, based on silhouette scores
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "figures/2_subtypes/kmean_k_silhouettes.pdf",
        "figures/2_subtypes/kmean_k_cosine.pdf",
        "figures/2_subtypes/kmean_k_sd.pdf",
        "figures/2_subtypes/kmean_k_abs_cosine.pdf"

    log:
        "logs/evaluate_kmeans_k.log"

    shell:
        "Rscript bin/analysis/2_subtypes/evaluate_kmeans_k.R &> {log}"

rule compare_hclust_dynamic_deep_split:
    """
    Assess different values of the cutreeHybrid deepSplit parameter
    """
    input:
        "data/combined_mutational_scans.tsv",
        "data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0_no_permissive.tsv",
        "data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1_no_permissive.tsv"

    output:
        "figures/2_subtypes/hclust_dynamic_deep_split_mean_corelation.pdf"

    log:
        "logs/compare_hclust_dynamic_deep_split.log"

    shell:
        "Rscript bin/analysis/2_subtypes/compare_hclust_dynamic_deep_split.R &> {log}"


rule compare_subtypes:
    """
    Compare clustering performance across the tested methods
    """
    input:
        subtypes=[f"data/subtypes/{c}.tsv" for c in STANDARD_CLUSTERINGS],
        dms="data/combined_mutational_scans.tsv"

    output:
        "figures/2_subtypes/method_silhouettes.pdf",
        "figures/2_subtypes/method_cosine_sim.pdf"

    log:
        "logs/compare_subtypes.log"

    shell:
        "Rscript bin/analysis/2_subtypes/compare_subtypes.R --dms {input.dms} --figures figures/2_subtypes --subtypes {input.subtypes} &> {log}"

rule final_subtypes:
    """
    Additional analysis of the final chosen subtypes
    """
    input:
        sections=ancient("meta/structures.yaml"),
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_0.tsv",
        alt_subtypes="data/subtypes/hclust_pca_no_sig_dynamic_cos_deep_1.tsv",

    output:
        "figures/2_subtypes/final_subtypes/outlier_profiles.pdf",
        directory("figures/2_subtypes/final_subtypes/cluster_heatmaps")

    log:
        "logs/final_subtypes.log"

    shell:
        "Rscript bin/analysis/2_subtypes/final_subtypes.R &> {log}"


#### Clustering and Analysing aach cluster method ####
## Kmeans
rule make_subtypes_kmeans:
    """
    Make Kmeans subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/kmeans_{name}.yaml'

    output:
        'data/subtypes/kmeans_{name}.tsv',
        'data/subtypes/kmeans_{name}.rds',
        [f'figures/2_subtypes/kmeans_{{name}}/{x}' for x in diagnostic_plots]

    log:
        'logs/make_subtypes/kmeans_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/kmeans_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/kmeans_{wildcards.name} --out 'data/subtypes/kmeans_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_kmeans:
    """
    Make characterisation plots for each amino acid from input Kmeans subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/kmeans_{name}.tsv",
        rds="data/subtypes/kmeans_{name}.rds"

    output:
        [f'figures/2_subtypes/kmeans_{{name}}/{x}' for x in characterisation_plots]

    log:
        'logs/characterise_subtypes/kmeans_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/kmeans_{wildcards.name} data/subtypes/kmeans_{wildcards.name} &> {log}"

## PAM
rule make_subtypes_pam:
    """
    Make PAM subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/pam_{name}.yaml'

    output:
        'data/subtypes/pam_{name}.tsv',
        'data/subtypes/pam_{name}.rds',
        [f'figures/2_subtypes/pam_{{name}}/{x}' for x in diagnostic_plots],
        'figures/2_subtypes/pam_{name}/medoids.pdf'

    log:
        'logs/make_subtypes/pam_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/pam_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/pam_{wildcards.name} --out 'data/subtypes/pam_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_pam:
    """
    Make characterisation plots for each amino acid from input PAM subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/pam_{name}.tsv",
        rds="data/subtypes/pam_{name}.rds"

    output:
        [f'figures/2_subtypes/pam_{{name}}/{x}' for x in characterisation_plots]

    log:
        'logs/characterise_subtypes/pam_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/pam_{wildcards.name} data/subtypes/pam_{wildcards.name} &> {log}"

## Hclust
rule make_subtypes_hclust:
    """
    Make Hclust subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/hclust_{name}.yaml'

    output:
        'data/subtypes/hclust_{name}.tsv',
        'data/subtypes/hclust_{name}.rds',
        [f'figures/2_subtypes/hclust_{{name}}/{x}' for x in diagnostic_plots],
        'figures/2_subtypes/hclust_{name}/hclust_dend.pdf'

    log:
        'logs/make_subtypes/hclust_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/hclust_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/hclust_{wildcards.name} --out 'data/subtypes/hclust_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_hclust:
    """
    Make characterisation plots for each amino acid from input Hclust subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/hclust_{name}.tsv",
        rds="data/subtypes/hclust_{name}.rds"

    output:
        [f'figures/2_subtypes/hclust_{{name}}/{x}' for x in characterisation_plots],
        [f'figures/2_subtypes/hclust_{{name}}/minimal_dends/{x}.pdf' for x in AA_ALPHABET]

    log:
        'logs/characterise_subtypes/hclust_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/hclust_{wildcards.name} data/subtypes/hclust_{wildcards.name} &> {log}"

## HDBSCAN
rule make_subtypes_hdbscan:
    """
    Make HDBSCAN subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/hdbscan_{name}.yaml'

    output:
        'data/subtypes/hdbscan_{name}.tsv',
        'data/subtypes/hdbscan_{name}.rds',
        [f'figures/2_subtypes/hdbscan_{{name}}/{x}' for x in diagnostic_plots],
        'figures/2_subtypes/hdbscan_{name}/hdbscan_dend.pdf'

    log:
        'logs/make_subtypes/hdbscan_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/hdbscan_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/hdbscan_{wildcards.name} --out 'data/subtypes/hdbscan_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_hdbscan:
    """
    Make characterisation plots for each amino acid from input HDBSCAN subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/hdbscan_{name}.tsv",
        rds="data/subtypes/hdbscan_{name}.rds"

    output:
        [f'figures/2_subtypes/hdbscan_{{name}}/{x}' for x in characterisation_plots]

    log:
        'logs/characterise_subtypes/hdbscan_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/hdbscan_{wildcards.name} data/subtypes/hdbscan_{wildcards.name} &> {log}"

## DBSCAN
rule make_subtypes_dbscan:
    """
    Make DBSCAN subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/dbscan_{name}.yaml'

    output:
        'data/subtypes/dbscan_{name}.tsv',
        'data/subtypes/dbscan_{name}.rds',
        [f'figures/2_subtypes/dbscan_{{name}}/{x}' for x in diagnostic_plots]

    log:
        'logs/make_subtypes/dbscan_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/dbscan_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/dbscan_{wildcards.name} --out 'data/subtypes/dbscan_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_dbscan:
    """
    Make characterisation plots for each amino acid from input DBSCAN subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/dbscan_{name}.tsv",
        rds="data/subtypes/dbscan_{name}.rds"

    output:
        [f'figures/2_subtypes/dbscan_{{name}}/{x}' for x in characterisation_plots]

    log:
        'logs/characterise_subtypes/dbscan_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/dbscan_{wildcards.name} data/subtypes/dbscan_{wildcards.name} &> {log}"

## GMM
rule make_subtypes_gmm:
    """
    Make GMM subtypes, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/gmm_{name}.yaml'

    output:
        'data/subtypes/gmm_{name}.tsv',
        'data/subtypes/gmm_{name}.rds',
        [f'figures/2_subtypes/gmm_{{name}}/{x}' for x in diagnostic_plots],
        'figures/2_subtypes/gmm_{name}/gmm_bic.pdf'

    log:
        'logs/make_subtypes/gmm_{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/gmm_{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/gmm_{wildcards.name} --out 'data/subtypes/gmm_{wildcards.name}' {input.yaml} &> {log}
        """

rule characterise_subtypes_gmm:
    """
    Make characterisation plots for each amino acid from input GMM subtypes
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/gmm_{name}.tsv",
        rds="data/subtypes/gmm_{name}.rds"

    output:
        [f'figures/2_subtypes/gmm_{{name}}/{x}' for x in characterisation_plots]

    log:
        'logs/characterise_subtypes/gmm_{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/gmm_{wildcards.name} data/subtypes/gmm_{wildcards.name} &> {log}"