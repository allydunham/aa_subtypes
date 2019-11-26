"""
Rules for SIFT pipeline

Expects global variables:
- GENES: dict mapping genes to lists of the studies assessing them
"""
FASTA_LINE_LENGTH = 80

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
    """
    Run SIFT4G on a FASTA file, assessing all possible variants.
    Note: for this project I have been using a modified version of SIFT4G that
    outputs to 4.d.p rather than 2.
    """
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
