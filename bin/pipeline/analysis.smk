"""
Rules for analysing combined data and calculating AA subtypes

Expects global variables:
- STANDARD_CLUSTERINGS: list of names for clusterings generated
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
        'figures/1_landscape/pc1_vs_mean_sift.pdf',
        'figures/1_landscape/pc1_pc2_mean_score.pdf',
        'figures/1_landscape/pc2_pc4_surface_accessibility.pdf',
        'figures/1_landscape/pc2_vs_surface_accessibility.pdf',
        'figures/1_landscape/pc2_vs_hydrophobicity.pdf',
        'figures/1_landscape/pc3_foldx_terms.pdf',
        'figures/1_landscape/foldx_pc_cor.pdf',
        'figures/1_landscape/tsne_study.pdf',
        'figures/1_landscape/tsne_aa.pdf',
        'figures/1_landscape/tsne_hydrophobicity.pdf',
        'figures/1_landscape/tsne_mean_er.pdf',
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

rule project_landscape:
    """
    Project landscape factors onto protein structures
    """
    input:
        dms='data/combined_mutational_scans.tsv',
        pdb='data/pdb/{gene}.pdb',
        yaml=ancient('meta/structures.yaml')

    output:
        [f'figures/1_landscape/pdb/{{gene}}/{{gene}}_{property}.png' for property in PDB_LANDSCAPE_FACTORS]

    log:
        'logs/project_landscape/{gene}.log'

    shell:
        f"python bin/analysis/1_landscape/project_landscape.py --pdb {{input.pdb}} --gene {{wildcards.gene}} --output_dir figures/1_landscape/pdb/{{wildcards.gene}} --structure_yaml {{input.yaml}} --data {{input.dms}} {' '.join(PDB_LANDSCAPE_FACTORS)} &> {{log}}"

rule project_landscape_colourbar:
    """
    Project landscape factors onto protein structures
    """
    input:
        dms='data/combined_mutational_scans.tsv'

    output:
        'figures/1_landscape/pdb/{property}_colourbar.pdf'

    log:
        'logs/project_landscape_colourbar/{property}.log'

    group:
        'project_colourbars'

    shell:
        'python bin/analysis/1_landscape/project_landscape.py --colourbar --output_dir figures/1_landscape/pdb --data {input.dms} {wildcards.property} &> {log}'

#### Clustering ####
diagnostic_plots = ['clustering.pdf', 'umap.pdf', 'tsne.pdf',
                    'silhouette_global.pdf', 'silhouette_clustering_vars.pdf',
                    'silhouette_per_aa.pdf', 'silhouette_per_aa_cosine.pdf']
rule make_subtypes:
    """
    Generalised clustering script, taking parameters from YAML files in meta/clustering
    """
    input:
        dms_wide='data/combined_mutational_scans.tsv',
        yaml='meta/subtypes/{name}.yaml'

    output:
        'data/subtypes/{name}.tsv',
        'data/subtypes/{name}.rds',
        [f'figures/2_subtypes/{{name}}/{x}' for x in diagnostic_plots]

    log:
        'logs/make_subtypes/{name}.log'

    shell:
        """
        mkdir figures/2_subtypes/{wildcards.name} &> {log} || true
        Rscript bin/analysis/2_subtypes/make_subtypes.R --figures figures/2_subtypes/{wildcards.name} {input.yaml} --out 'data/subtypes/{wildcards.name}' &> {log}
        """

characterisation_plots = ['ramachandran_angles.pdf', 'sizes.pdf', 'er_profiles.pdf',
                          'er_correlation.pdf', 'foldx.pdf', 'chem_env.pdf', 'ss_probability.pdf',
                          'er_vs_surface_accessibility.pdf', 'er_vs_size.pdf']
rule characterise_subtypes:
    """
    Make characterisation plots for each amino acid from an input subtype clustering
    """
    input:
        dms="data/combined_mutational_scans.tsv",
        subtypes="data/subtypes/{name}.tsv"

    output:
        [f'figures/2_subtypes/{{name}}/{x}' for x in characterisation_plots],
        [f"figures/2_subtypes/{{name}}/aa_profiles/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_subtypes/{{name}}/aa_profiles_relative/{x}.pdf" for x in AA_ALPHABET],
        [f"figures/2_subtypes/{{name}}/ss_probabilities/{x}.pdf" for x in AA_ALPHABET]

    log:
        'logs/characterise_subtypes/{name}.log'

    shell:
        "Rscript bin/analysis/2_subtypes/characterise_subtypes.R --dms {input.dms} --figures figures/2_subtypes/{wildcards.name} {input.subtypes} &> {log}"

rule all_position_subtypes:
    """
    Genreate subtypes from positions of all amino acids (using Hclust with dynamic cutting)
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "data/clustering/hclust_profile_dynamic_all_positions.tsv",
        [f"figures/2_subtypes/hclust_profile_dynamic_all_positions/{x}" for x in diagnostic_plots],
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