"""
Rules for analysing combined mutational landscape
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

rule continuous_characterisation:
    """
    Analyse all positions against the various biochemical factors directly
    """
    input:
        "data/combined_mutational_scans.tsv"

    output:
        "figures/3_continuous/mean_er_correlations.pdf"

    log:
        "logs/continuous_characterisation.log"

    shell:
        "Rscript bin/analysis/3_continuous/all_position_characterisation.R &> {log}"