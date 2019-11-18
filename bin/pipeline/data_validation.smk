"""
Rules for validating study data

Expects global variables:
- STUDIES: dictionary ofall study dicts (as read from the yamls)
- GENES: dictionary mapping genes to lists of their study IDs
"""

VALIDATION_PLOTS = [
    "figures/0_data_properties/sift_score_correlation.pdf",
    "figures/0_data_properties/per_study/melnikov_2014_aph3ii/initial_library_correlation.pdf",
    "figures/0_data_properties/per_study/melnikov_2014_aph3ii/filtered_library_correlation.pdf",
    "figures/0_data_properties/per_study/melnikov_2014_aph3ii/rel_conc_correlation.pdf",
    "figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_correlation.pdf",
    "figures/0_data_properties/per_study/kitzman_2015_gal4/validate_selection_combination.pdf",
    "figures/0_data_properties/per_study/giacomelli_2018_tp53/initial_experiment_cor.pdf",
    "figures/0_data_properties/per_study/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf",
    "figures/0_data_properties/per_study/giacomelli_2018_tp53/conditions.pdf",
    "figures/0_data_properties/per_study/heredia_2018_ccr5/replicate_correlation.pdf",
    "figures/0_data_properties/per_study/heredia_2018_ccr5/experiment_correlation.pdf",
    "figures/0_data_properties/per_study/heredia_2018_cxcr4/replicate_correlation.pdf",
    "figures/0_data_properties/per_study/heredia_2018_cxcr4/experiment_correlation.pdf",
    "figures/0_data_properties/per_study/sarkisyan_2016_gfp/multi_mut_validation.pdf",
    "figures/0_data_properties/per_study/dorrity_2018_ste12/rep_correlation.pdf",
    "figures/0_data_properties/per_study/dorrity_2018_ste12/multi_mut_validation.pdf",
    "figures/0_data_properties/per_study/araya_2012_yap1/multi_mut_validation.pdf",
    "figures/0_data_properties/per_study/starita_2013_ube4b/multi_mut_validation.pdf",
    "figures/0_data_properties/sift_score_correlation.pdf",
    "figures/0_data_properties/sift_score_density.pdf",
    "figures/0_data_properties/gene_repeats/brca1.pdf",
    "figures/0_data_properties/gene_repeats/hsp90.pdf",
    "figures/0_data_properties/gene_repeats/tem1.pdf",
    "figures/0_data_properties/gene_repeats/ubi.pdf"]

rule validate_melnikov:
    input:
        [f'data/studies/melnikov_2014_aph3ii/raw/{x}' for
         x in STUDIES['melnikov_2014_aph3ii']['input_files']]

    output:
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/initial_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/filtered_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/rel_conc_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_correlation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_melnikov_2014_aph3ii.R"

# Validate Kitzman et al. 2015 (GAL4)
rule validate_kitzman:
    input:
        [f'data/studies/kitzman_2015_gal4/raw/{x}' for
         x in STUDIES['kitzman_2015_gal4']['input_files']]

    output:
        "figures/0_data_properties/per_study/kitzman_2015_gal4/validate_selection_combination.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_kitzman_2015_gal4.R"

# Validate Giacomelli et al. 2018 (TP53)
rule validate_giacomelli:
    input:
        [f'data/studies/giacomelli_2018_tp53/raw/{x}' for
         x in STUDIES['giacomelli_2018_tp53']['input_files']]

    output:
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/initial_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/conditions.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_giacomelli_2018_tp53.R"

# Validate Heredia et al. 2018
rule validate_heredia:
    input:
        [f'data/studies/heredia_2018_ccr5/raw/{x}' for
         x in STUDIES['heredia_2018_ccr5']['input_files']] +
        [f'data/studies/heredia_2018_cxcr4/raw/{x}' for
         x in STUDIES['heredia_2018_cxcr4']['input_files']]

    output:
        "figures/0_data_properties/per_study/heredia_2018_ccr5/replicate_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_ccr5/experiment_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_cxcr4/replicate_correlation.pdf",
        "figures/0_data_properties/per_study/heredia_2018_cxcr4/experiment_correlation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_heredia_2018.R"

# Validate Sarkisyan et al. 2016 (GFP)
rule validate_sarkisyan:
    input:
        [f'data/studies/sarkisyan_2016_gfp/raw/{x}' for
         x in STUDIES['sarkisyan_2016_gfp']['input_files']]

    output:
        "figures/0_data_properties/per_study/sarkisyan_2016_gfp/multi_mut_validation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_sarkisyan_2016_gfp.R"

# Validate Dorrity et al. 2018 (STE12)
rule validate_dorrity:
    input:
        [f'data/studies/dorrity_2018_ste12/raw/{x}' for
         x in STUDIES['dorrity_2018_ste12']['input_files']]

    output:
        "figures/0_data_properties/per_study/dorrity_2018_ste12/rep_correlation.pdf",
        "figures/0_data_properties/per_study/dorrity_2018_ste12/multi_mut_validation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_dorrity_2018_ste12.R"

# Validate Araya et al. 2012 (YAP1)
rule validate_araya:
    input:
        [f'data/studies/araya_2012_yap1/raw/{x}' for
         x in STUDIES['araya_2012_yap1']['input_files']]

    output:
        "figures/0_data_properties/per_study/araya_2012_yap1/multi_mut_validation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_araya_2012_yap1.R"

# Validate Starita et al. 2013 (UBE4B)
rule validate_starita:
    input:
        [f'data/studies/starita_2013_ube4b/raw/{x}' for
         x in STUDIES['starita_2013_ube4b']['input_files']]

    output:
        "figures/0_data_properties/per_study/starita_2013_ube4b/multi_mut_validation.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_starita_2013_ube4b.R"

# Check correlation with SIFT
rule sift_correlation:
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv')),
        expand('data/sift/{gene}.{ext}', gene=GENES.keys(), ext=('fa', 'SIFTprediction'))

    output:
        'figures/0_data_properties/sift_score_correlation.pdf',
        'figures/0_data_properties/sift_score_density.pdf'

    shell:
        "Rscript bin/analysis/0_data_properties/sift_correlation.R"

# Analyse cases with multiple studies of the same gene
rule gene_repeats:
    input:
        expand('data/studies/{study}/{study}.{ext}',
               ext=['yaml', 'tsv'],
               study=['findlay_2018_brca1', 'starita_2015_brca1',
                      'hietpas_2011_hsp90', 'jiang_2013_hsp90',
                      'mishra_2016_hsp90', 'firnberg_2014_tem1',
                      'steinberg_2016_tem1', 'roscoe_2013_ubi',
                      'roscoe_2014_ubi'])

    output:
        "figures/0_data_properties/gene_repeats/brca1.pdf",
        "figures/0_data_properties/gene_repeats/hsp90.pdf",
        "figures/0_data_properties/gene_repeats/tem1.pdf",
        "figures/0_data_properties/gene_repeats/ubi.pdf"

    shell:
        "Rscript bin/analysis/0_data_properties/gene_repeats.R"

