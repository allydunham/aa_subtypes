"""
Rules for generating figures
"""

rule figure1:
    """
    Figure 1
    """
    input:
        "data/combined_mutational_scans.tsv",
        expand('data/studies/{study}/{study}.{ext}',
               study=UNFILTERED_STUDIES, ext=('yaml', 'tsv')),
        "data/studies/roscoe_2013_ubi/roscoe_2013_ubi.tsv",
        "data/studies/roscoe_2013_ubi/roscoe_2013_ubi.yaml",
        "data/studies/roscoe_2014_ubi/roscoe_2014_ubi.tsv",
        "data/studies/roscoe_2014_ubi/roscoe_2014_ubi.yaml"

    output:
        "figures/4_figures/figure1.pdf",
        "figures/4_figures/figure1.png"

    log:
        "logs/figure1.log"

    shell:
        "Rscript bin/figures/figure1.R &> {log}"

rule figure2:
    """
    Figure 2
    """
    input:
        "data/combined_mutational_scans.tsv",
        "meta/uniprot_domains.gff",
        [f"figures/4_figures/proteins/{p}.png" for p in UNFILTERED_GENES]

    output:
        "figures/4_figures/figure2.pdf",
        "figures/4_figures/figure2.png"

    log:
        "logs/figure2.log"

    shell:
        "Rscript bin/figures/figure2.R &> {log}"

rule figure3:
    """
    Figure 3
    """
    input:
        "data/combined_mutational_scans.tsv",
        "data/subtypes/final_subtypes.tsv",

    output:
        "figures/4_figures/figure3.pdf",
        "figures/4_figures/figure3.png",
        "figures/4_figures/parts/figure3_cor_heatmap.pdf",
        "figures/4_figures/parts/figure3_cor_dendrogram.pdf",
        "figures/4_figures/parts/figure3_cor_set_small_aliphatic.pdf",
        "figures/4_figures/parts/figure3_cor_set_not_proline.pdf",
        "figures/4_figures/parts/figure3_cor_set_positive.pdf",
        "figures/4_figures/parts/figure3_cor_set_aromatic.pdf",
        "figures/4_figures/parts/figure3_cor_set_aliphatic.pdf",
        "figures/4_figures/parts/figure3_cor_set_larger_aliphatic.pdf",
        "figures/4_figures/parts/figure3_cor_set_not_aromatic.pdf",
        "figures/4_figures/parts/figure3_cor_set_negative.pdf"

    log:
        "logs/figure3.log"

    shell:
        "Rscript bin/figures/figure3.R &> {log}"

rule figure4:
    """
    Figure 4
    """
    input:
        "data/combined_mutational_scans.tsv",
        "data/subtypes/final_subtypes.tsv"

    output:
        "figures/4_figures/figure4.pdf",
        "figures/4_figures/figure4.png"

    log:
        "logs/figure4.log"

    shell:
        "Rscript bin/figures/figure4.R &> {log}"