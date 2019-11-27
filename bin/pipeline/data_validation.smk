"""
Rules for validating study data

Expects global variables:
- STUDIES: dictionary ofall study dicts (as read from the yamls)
- GENES: dictionary mapping genes to lists of their study IDs
"""

#### Validate standardiation methods ####
rule validate_melnikov:
    """
    Validate standardisation method used for Melnikov at al. 2014 (APH(3')II)
    """
    input:
        [f'data/studies/melnikov_2014_aph3ii/raw/{x}' for
         x in STUDIES['melnikov_2014_aph3ii']['input_files']]

    output:
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/initial_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/filtered_library_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/rel_conc_correlation.pdf",
        "figures/0_data_properties/per_study/melnikov_2014_aph3ii/drug_correlation.pdf"

    log:
        "logs/data_validation/melnikov_2014_aph3ii.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_melnikov_2014_aph3ii.R &> {log}"

rule validate_kitzman:
    """
    Validate standardisation method used for Kitzman et al. 2015 (GAL4)
    """
    input:
        [f'data/studies/kitzman_2015_gal4/raw/{x}' for
         x in STUDIES['kitzman_2015_gal4']['input_files']]

    output:
        "figures/0_data_properties/per_study/kitzman_2015_gal4/validate_selection_combination.pdf"

    log:
        "logs/data_validation/kitzman_2015_gal4.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_kitzman_2015_gal4.R &> {log}"

rule validate_giacomelli:
    """
    Validate standardisation method used for Giacomelli et al. 2018 (TP53)
    """
    input:
        [f'data/studies/giacomelli_2018_tp53/raw/{x}' for
         x in STUDIES['giacomelli_2018_tp53']['input_files']]

    output:
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/initial_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/codon_averaged_experiment_cor.pdf",
        "figures/0_data_properties/per_study/giacomelli_2018_tp53/conditions.pdf"

    log:
        "logs/data_validation/giacomelli_2018_tp53.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_giacomelli_2018_tp53.R &> {log}"

rule validate_heredia:
    """
    Validate standardisation method used for Heredia et al. 2018 (CCR5 & CXCR4)
    """
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

    log:
        "logs/data_validation/heredia_2018.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_heredia_2018.R &> {log}"

rule validate_sarkisyan:
    """
    Validate standardisation method used for Sarkisyan et al. 2016 (GFP)
    """
    input:
        [f'data/studies/sarkisyan_2016_gfp/raw/{x}' for
         x in STUDIES['sarkisyan_2016_gfp']['input_files']]

    output:
        "figures/0_data_properties/per_study/sarkisyan_2016_gfp/multi_mut_validation.pdf"

    log:
        "logs/data_validation/sarkisyan_2016_gfp.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_sarkisyan_2016_gfp.R &> {log}"

rule validate_dorrity:
    """
    Validate standardisation method used for Dorrity et al. 2018 (STE12)
    """
    input:
        [f'data/studies/dorrity_2018_ste12/raw/{x}' for
         x in STUDIES['dorrity_2018_ste12']['input_files']]

    output:
        "figures/0_data_properties/per_study/dorrity_2018_ste12/rep_correlation.pdf",
        "figures/0_data_properties/per_study/dorrity_2018_ste12/multi_mut_validation.pdf"

    log:
        "logs/data_validation/dorrity_2018_ste12.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_dorrity_2018_ste12.R &> {log}"

rule validate_araya:
    """
    Validate standardisation method used for Araya et al. 2012 (YAP1)
    """
    input:
        [f'data/studies/araya_2012_yap1/raw/{x}' for
         x in STUDIES['araya_2012_yap1']['input_files']]

    output:
        "figures/0_data_properties/per_study/araya_2012_yap1/multi_mut_validation.pdf"

    log:
        "logs/data_validation/araya_2012_yap1.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_araya_2012_yap1.R &> {log}"

rule validate_starita:
    """
    Validate standardisation method used for Starita et al. 2013 (UBE4B)
    """
    input:
        [f'data/studies/starita_2013_ube4b/raw/{x}' for
         x in STUDIES['starita_2013_ube4b']['input_files']]

    output:
        "figures/0_data_properties/per_study/starita_2013_ube4b/multi_mut_validation.pdf"

    log:
        "logs/data_validation/starita_2013_ube4b.log"

    shell:
        "Rscript bin/analysis/0_data_properties/validate_starita_2013_ube4b.R &> {log}"

#### Validate overall dataset ####
rule sift_correlation:
    """
    Check correlation between study results and SIFT scores
    """
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv')),
        expand('data/sift/{gene}.{ext}', gene=GENES.keys(), ext=('fa', 'SIFTprediction'))

    output:
        'figures/0_data_properties/sift_score_correlation.pdf',
        'figures/0_data_properties/sift_score_density.pdf'

    log:
        "logs/data_validation/sift_correlation.log"

    shell:
        "Rscript bin/analysis/0_data_properties/sift_correlation.R &> {log}"

# Analyse cases with multiple studies of the same gene
rule gene_repeats:
    """
    Analyse correlation between multiple studies of the same genes
    """
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

    log:
        "logs/data_validation/gene_repeats.log"

    shell:
        "Rscript bin/analysis/0_data_properties/gene_repeats.R &> {log}"

VALIDATION_RULES = [rules.validate_melnikov, rules.validate_kitzman, rules.validate_giacomelli,
                    rules.validate_heredia, rules.validate_sarkisyan, rules.validate_dorrity,
                    rules.validate_araya, rules.validate_starita, rules.sift_correlation,
                    rules.gene_repeats]

rule validate_data:
    """
    Produce all data validation plots
    """
    input:
        [r.output for r in VALIDATION_RULES]
