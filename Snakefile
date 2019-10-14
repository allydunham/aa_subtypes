# pylint: disable-all
"""
Pipeline for the Mutational Landscapes/Amino Acids Subtypes Project
"""

from ruamel.yaml import YAML

yaml = YAML(typ='safe')

UNIREF90_DB_PATH = '/hps/research1/beltrao/ally/databases/uniref90/uniref90_2019_1.fasta'

# Process the raw data from each study
rule standardise_data:
    input:
        "data/studies/{study}/standardise_{study}.R"

    output:
        "data/studies/{study}/{study}.tsv"

    shell:
        "Rscript {input}"

# Make all SIFT predictions for study genes
# TODO Change layout so as to only run each gene once
rule make_fasta:
    input:
        "data/studies/{study}/{study}.yaml"

    output:
        "data/studies/{study}/{study}.fa"

    shell:
        "python bin/make_study_fasta.py -l 80 {input} > {output}"

rule sift4g:
    input:
        fa = "data/studies/{study}/{study}.fa",
        db = UNIREF90_DB_PATH

    output:
        "data/studies/{study}/{study}.SIFTprediction"

    shell:
        "sift4g -q {input.fa} -d {input.db} --out data/studies/{study}"

# Make all FoldX for study genes
# TODO