"""
Rules for generating figures
"""

figure_inputs = {
    '1': expand('data/studies/{study}/{study}.{ext}',
                study=UNFILTERED_STUDIES, ext=('yaml', 'tsv')) +
         ["data/combined_mutational_scans.tsv", "figures/4_figures/parts/arrow.png"],

    '2': [f"figures/4_figures/proteins/{p}.png" for p in UNFILTERED_GENES] +
         [ "data/combined_mutational_scans.tsv", "meta/uniprot_domains.gff"],

    '3': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"],

    '4': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"] +
         ["figures/4_figures/position_examples/gal4_cys_zinc.png",
          "figures/4_figures/position_examples/cbs_cys_haem.png",
          "figures/4_figures/position_examples/np_cys_aromatic.png",
          "figures/4_figures/position_examples/ccr5_cys_aromatic.png"],

    '5': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"] +
         ["figures/4_figures/position_examples/cbs_asp_ionic.png",
          "figures/4_figures/position_examples/tem1_asp_ligand.png",
          "figures/4_figures/position_examples/gal4_lys_dna.png",
          "figures/4_figures/position_examples/tem1_asp_sa.png",
          "figures/4_figures/position_examples/aph3ii_arg_not_proline.png",
          "figures/4_figures/position_examples/pab1_arg_not_neg.png"],

    '6': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"] +
         ["figures/4_figures/position_examples/ras_aliphatic_entropy.png",
          "figures/4_figures/position_examples/adrb2_ala_small_hydro.png",
          "figures/4_figures/position_examples/ras_met_buried.png",
          "figures/4_figures/position_examples/cbs_phe_pi.png",
          "figures/4_figures/position_examples/tem1_trp_surface.png"],

    'S1': ["data/studies/araya_2012_yap1/araya_2012_yap1.yaml",
           "data/studies/araya_2012_yap1/raw/urn_mavedb_00000002-a-2_scores.csv",
           "data/studies/dorrity_2018_ste12/dorrity_2018_ste12.yaml",
           "data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd01.xlsx",
           "data/studies/starita_2013_ube4b/starita_2013_ube4b.yaml",
           "data/studies/starita_2013_ube4b/raw/sd01.xlsx"],

    'S2': expand('data/studies/{study}/{study}.{ext}', ext=['yaml', 'tsv'],
                 study=['findlay_2018_brca1', 'starita_2015_brca1',
                        'hietpas_2011_hsp90', 'jiang_2013_hsp90',
                        'firnberg_2014_tem1', 'steinberg_2016_tem1',
                        'roscoe_2013_ubi', 'roscoe_2014_ubi']),

    'S3': "data/combined_mutational_scans.tsv",

    'S4': "data/combined_mutational_scans.tsv",

    'S5': "data/combined_mutational_scans.tsv",

    'S6': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"],

    'S7': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"],

    'S28': ["meta/uniprot_domains.gff", "data/long_combined_mutational_scans.tsv",
            "data/combined_mutational_scans.tsv"],

    'S29': ["data/subtypes/final_subtypes.tsv", "data/subtypes/sift_scores.tsv",
            "data/combined_mutational_scans.tsv"]
}

rule figureS8_27:
    """
    Generate subtype characterisation figures for each amino acids subtypes (figures S8-27)
    """
    input:
        "data/combined_mutational_scans.tsv",
        "data/subtypes/final_subtypes.tsv"

    output:
        expand("figures/4_figures/figureS{i}.{ext}", i=range(8,28), ext=['pdf', 'png']),
        "figures/4_figures/figureS8_27.pdf"

    log:
        "logs/figureS8_27.log"

    shell:
        "Rscript bin/figures/figureS8_27.R &> {log}"

rule figure:
    """
    Generate a numbered figure for the paper
    """
    input: lambda w: figure_inputs[w.num]

    output:
        "figures/4_figures/figure{num}.pdf",
        "figures/4_figures/figure{num}.png"

    wildcard_constraints:
        num="1|2|3|4|5|6|S1|S2|S3|S4|S5|S6|S7|S28|S29"

    log:
        "logs/figure{num}.log"

    shell:
        "Rscript bin/figures/figure{wildcards.num}.R &> {log}"


