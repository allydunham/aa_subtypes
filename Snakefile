# pylint: disable-all
"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""
import os
import math
from collections import defaultdict
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from ruamel.yaml import YAML
import subtypes_utils as sutil

configfile: 'snakemake.yaml'
localrules: all, clean, all_standardisation, all_sift, make_gene_fasta,
    foldx_variants, foldx_split, foldx_combine

yaml = YAML(typ='safe')

AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'

# Hash of study IDs to their config
STUDIES = {}
for study in os.listdir('data/studies'):
    with open(f'data/studies/{study}/{study}.yaml') as yaml_file:
        STUDIES[study] = yaml.load(yaml_file)

# Hash linking genes to studies
GENES = defaultdict(list)
for study, conf in STUDIES.items():
    GENES[sutil.gene_to_filename(conf['gene'])].append(study)

# Explicitly note all the output here
# even though some would necessarily cause other bits to generate
rule all:
    input:
        expand('data/studies/{study}/{study}.tsv', study=STUDIES.keys()),
        expand('data/sift/{gene}.fa', gene=GENES.keys()),
        expand('data/sift/{gene}.SIFTprediction', gene=GENES.keys()),
        expand("data/foldx/{gene}/{gene}_Repair.pdb", gene=GENES.keys()),
        expand("data/foldx/{gene}/individual_list", gene=GENES.keys()),
        expand("data/foldx/{gene}/average_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/dif_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/raw_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/pdblist_{gene}.fxout", gene=GENES.keys()),
        'meta/study_summary.tsv',
        'meta/gene_summary.tsv',
        'meta/overall_summary',
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
        'figures/0_data_properties/sift_score_correlation.pdf',
        'figures/0_data_properties/sift_score_density.pdf',
        'figures/0_data_properties/study_variants_summary.pdf',
        'figures/0_data_properties/gene_variants_summary.pdf',
        "figures/0_data_properties/gene_repeats/brca1.pdf",
        "figures/0_data_properties/gene_repeats/hsp90.pdf",
        "figures/0_data_properties/gene_repeats/tem1.pdf",
        "figures/0_data_properties/gene_repeats/ubi.pdf",
        'figures/0_data_properties/position_coverage.pdf'


# TODO Protection for files that take a long time to build? (Mainly SIFT for now)
# TODO more specific cleaning of sift dir
rule clean:
    run:
        output_files = [f'data/studies/{s}/{s}.tsv' for s in STUDIES.keys()]
        output_files.append('-r figures/0_data_properties/*')
        output_files.append('data/sift/*')
        output_files.append('meta/study_summary.tsv')
        output_files.append('meta/gene_summary.tsv')
        output_files.append('meta/overall_summary')

        for i in output_files:
            shell(f'rm {i} && echo "rm {i}" || true')

rule all_standardisation:
    input:
        expand('data/studies/{study}/{study}.tsv', study=STUDIES.keys())

rule all_sift:
    input:
        expand('data/sift/{gene}.SIFTprediction', gene=GENES.keys())

rule all_foldx:
    input:
        expand("data/foldx/{gene}/average_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/dif_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/raw_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/pdblist_{gene}.fxout", gene=GENES.keys())

#### Validate Data ####
# Validate Melnikov et al. 2014 (APH(3')-II)
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

# Summary plots of coverage/completeness
rule data_summary_plots:
    input:
        'meta/study_summary.tsv',
        'meta/gene_summary.tsv',
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        'figures/0_data_properties/study_variants_summary.pdf',
        'figures/0_data_properties/gene_variants_summary.pdf',
        'figures/0_data_properties/position_coverage.pdf'

    shell:
        'Rscript bin/analysis/0_data_properties/data_summary_plots.R'

#### Standardise Data ####
# Process the raw data from each study
rule standardise_study:
    input:
        "data/studies/{study}/standardise_{study}.R",
        lambda wildcards: [f'data/studies/{wildcards.study}/raw/{x}' for x in
                           STUDIES[wildcards.study]['input_files']]

    output:
        "data/studies/{study}/{study}.tsv",
        "figures/0_data_properties/per_study/{study}/original_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/transformed_distribution.pdf",
        "figures/0_data_properties/per_study/{study}/normalised_distribution.pdf"

    log:
        "logs/standardise_study/{study}.log"

    shell:
        "Rscript {input} 2> {log}"

rule summarise_studies:
    input:
        expand('data/studies/{study}/{study}.{ext}', study=STUDIES.keys(), ext=('yaml', 'tsv'))

    output:
        study='meta/study_summary.tsv',
        gene='meta/gene_summary.tsv',
        overall='meta/overall_summary',

    shell:
        "python bin/utils/summarise_studies.py -s {output.study} -g {output.gene} -u {output.overall} data/studies/*"

#### SIFT4G Predictions ####
rule make_gene_fasta:
    # seq really shouldn't change much even when .yaml's do -> sift results wont either
    # force re-run if neccessary by deleting relevant .fa
    input:
        lambda wildcards: [ancient(f'data/studies/{i}/{i}.yaml') for i
                           in GENES[wildcards.gene]]

    output:
        "data/sift/{gene}.fa"

    run:
        seq = None
        for study_yaml in input:
            with open(study_yaml, 'r') as yaml_file:
                conf = yaml.load(yaml_file)

            if seq is None:
                seq = conf['seq']
            elif not seq == conf['seq']:
                raise ValueError(f"Two studies on {wildcards.gene} have different sequences")

        with open(output[0], 'w') as fasta_file:
            print(f">{wildcards.gene}", file=fasta_file)
            for i in range(0, len(seq), FASTA_LINE_LENGTH):
                print(seq[i:(i + FASTA_LINE_LENGTH)], file=fasta_file)

rule sift4g:
    input:
        fa = "data/sift/{gene}.fa",
        db = config['sift']['uniref90_fa_path']

    output:
        "data/sift/{gene}.SIFTprediction"

    log:
        'logs/sift4g/{gene}.log'

    resources:
        mem_mb = 8000

    shell:
        "sift4g -q {input.fa} -d {input.db} --out data/sift 2> {log}"

#### FoldX5 Predictions ####
rule foldx_variants:
    input:
        pdb="data/foldx/{gene}/{gene}.pdb"

    output:
        muts="data/foldx/{gene}/individual_list"

    run:
        parser = PDBParser()
        structure = parser.get_structure(f'{wildcards.gene}', input.pdb)
        chain = structure[0]['A']

        variants = []
        for residue in chain:
            position = residue.id[1]
            aa = seq1(residue.get_resname())
            variants.extend([f'{aa}A{position}{x}' for x in AA_ALPHABET if not x == aa])

        with open(output.muts, 'w') as outfile:
            print(*variants, sep=';\n', end=';\n', file=outfile)

checkpoint foldx_split:
    input:
        "data/foldx/{gene}/individual_list"

    output:
        directory("data/foldx/{gene}/processing")

    params:
        n_lines = config['foldx']['variants_per_run']

    shell:
        """
        mkdir "data/foldx/{wildcards.gene}/processing"
        split -l {params.n_lines} data/foldx/{wildcards.gene}/individual_list data/foldx/{wildcards.gene}/processing/individual_list_
        """

rule foldx_repair:
    input:
        pdb="data/foldx/{gene}/{gene}.pdb"

    output:
        "data/foldx/{gene}/{gene}_Repair.pdb",
        "data/foldx/{gene}/{gene}_Repair.fxout"

    shell:
        "foldx --command=RepairPDB --pdb={input.pdb} --clean-mode=3"

rule foldx_model:
    input:
        pdb="data/foldx/{gene}/{gene}_Repair.pdb",
        muts="data/foldx/{gene}/processing/individual_list_{n}"

    output:
        "data/foldx/{gene}/processing/Average_{gene}_{n}_BM.fxout",
        "data/foldx/{gene}/processing/Dif_{gene}_{n}_BM.fxout",
        "data/foldx/{gene}/processing/Raw_{gene}_{n}_BM.fxout",
        "data/foldx/{gene}/processing/PdbList_{gene}_{n}_BM.fxout"

    shell:
        'FoldX --command=BuildModel --pdb={input.pdb} --mutant-file={input.muts} --output-file="{wildcards.gene}_{wildcards.n}" --output-dir=data/foldx/{wildcards.gene}/processing --numberOfRuns=3 --clean-mode=3 --out-pdb=false'

def get_foldx_split_files(wildcards):
    """
    Retrieve the split FoldX jobs
    """
    checkpoint_outdir = checkpoints.foldx_split.get(gene=wildcards.gene).output[0]
    return expand('data/foldx/{wildcards.gene}/processing/{fi}_{gene}_{n}_BM.fxout',
                   n=glob_wildcards(os.path.join(checkpoint_outdir, "individual_list_{n}")).n,
                   fi=('Average', 'Dif', 'Raw', 'PdbList'))

rule foldx_combine:
    input:
        get_foldx_split_files

    output:
        touch("data/foldx/{gene}/average_{gene}.fxout"),
        touch("data/foldx/{gene}/dif_{gene}.fxout"),
        touch("data/foldx/{gene}/raw_{gene}.fxout"),
        touch("data/foldx/{gene}/pdblist_{gene}.fxout")

    run:
        print('Add method for FoldX combining!')
