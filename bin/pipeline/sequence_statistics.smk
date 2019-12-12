"""
Rules for generating statistics from protein sequences, including the Porter5 and SIFT pipeline

Expects global variables:
- GENES: dict mapping genes to lists of the studies assessing them
"""

rule make_gene_fasta:
    """
    Generate a FASTA file for a gene from the sequence in the studies YAML configs
    """
    # marked ancient as seq shouldn't change when .yaml's do
    # force re-run if neccessary by deleting relevant .fa
    input:
        lambda wildcards: [ancient(f'data/studies/{i}/{i}.yaml') for i
                           in GENES[wildcards.gene]]

    output:
        "data/fasta/{gene}.fa"

    log:
        "logs/make_gene_fasta/{gene}.log"

    shell:
        "python bin/data_processing/make_gene_fasta.py {input} > {output} 2> {log}"

rule sift4g:
    """
    Run SIFT4G on a FASTA file, assessing all possible variants.
    Note: for this project I have been using a modified version of SIFT4G that
    outputs to 4.d.p rather than 2.
    """
    input:
        fa = "data/fasta/{gene}.fa",
        db = config['sift']['uniref90_fa_path']

    output:
        "data/sift/{gene}.SIFTprediction"

    log:
        'logs/sift4g/{gene}.log'

    resources:
        mem_mb = 8000

    shell:
        "sift4g -q {input.fa} -d {input.db} --out data/sift 2> {log}"

rule all_sift_predictions:
    """
    Produce SIFT predictions for all genes
    """
    input:
        expand('data/sift/{gene}.SIFTprediction', gene=GENES.keys())

rule porter5:
    """
    Run Porter5 on a FASTA file, generating secondary structure predictions for each residue.
    """
    input:
        "data/fasta/{gene}.fa"

    output:
        ss3="data/porter5/{gene}.ss3",
        ss8="data/porter5/{gene}.ss8"

    log:
        "logs/porter5/{gene}.log"

    threads: 4

    resources:
        mem_mb = 20000

    shell:
        f"""
        python {config['porter5']['path']} -i {{input}} --fast --cpu 4 --tmp &> {{log}}
        cat data/fasta/{{wildcards.gene}}.fa.log >> {{log}} 2>&1
        mv data/fasta/{{wildcards.gene}}.fa.ss3 {{output.ss3}} >> {{log}} 2>&1
        mv data/fasta/{{wildcards.gene}}.fa.ss8 {{output.ss8}} >> {{log}} 2>&1
        rm data/fasta/{{wildcards.gene}}.fa.* &> {{log}}
        rm data/fasta/{{wildcards.gene}}.hhr &> {{log}}
        """

rule all_porter5_predictions:
    """
    Produce SIFT predictions for all genes
    """
    input:
        expand('data/porter5/{gene}.ss8', gene=GENES.keys()),
        expand('data/porter5/{gene}.ss3', gene=GENES.keys())
