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
           "figures/4_figures/position_examples/cbs_cys_haem.png"
           "figures/4_figures/position_examples/np_cys_aromatic.png"
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

    'S1': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"],

    'S2': ["data/combined_mutational_scans.tsv", "data/subtypes/final_subtypes.tsv"],

    'S3': expand('data/studies/{study}/{study}.{ext}', ext=['yaml', 'tsv'],
                 study=['findlay_2018_brca1', 'starita_2015_brca1',
                        'hietpas_2011_hsp90', 'jiang_2013_hsp90',
                        'firnberg_2014_tem1', 'steinberg_2016_tem1',
                        'roscoe_2013_ubi', 'roscoe_2014_ubi']),

    'S4': "data/combined_mutational_scans.tsv",

    'S5': ["data/studies/araya_2012_yap1/araya_2012_yap1.yaml",
          "data/studies/araya_2012_yap1/raw/urn_mavedb_00000002-a-2_scores.csv",
          "data/studies/dorrity_2018_ste12/dorrity_2018_ste12.yaml",
          "data/studies/dorrity_2018_ste12/raw/pnas.1805882115.sd01.xlsx",
          "data/studies/starita_2013_ube4b/starita_2013_ube4b.yaml",
          "data/studies/starita_2013_ube4b/raw/sd01.xlsx"]
}

rule figure:
    """
    Generate a numbered figure for the paper
    """
    input: lambda w: figure_inputs[w.num]

    output:
        "figures/4_figures/figure{num}.pdf",
        "figures/4_figures/figure{num}.png"

    log:
        "logs/figure{num}.log"

    shell:
        "Rscript bin/figures/figure{wildcards.num}.R &> {log}"
