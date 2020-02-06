"""
Rules for FoldX pipeline

Expects one global variables:
- GENES: dict of all genes linking to lists of studies on each gene
"""

rule all_foldx_predictions:
    """
    Produce FoldX predictions for all genes
    """
    input:
        expand("data/foldx/{gene}/average_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/dif_{gene}.fxout", gene=GENES.keys()),
        expand("data/foldx/{gene}/raw_{gene}.fxout", gene=GENES.keys())

rule foldx_variants:
    """
    Produce a list of all possible variants from a PDB structure, from positions that
    correspond to the regions defined in meta/structures.yaml as part of the analyses protein
    """
    input:
        pdb="data/foldx/{gene}/{gene}_Repair.pdb",
        structures=ancient("meta/structures.yaml") # Imported generally

    output:
        muts="data/foldx/{gene}/individual_list"

    log:
        "logs/foldx_variants/{gene}.log"

    shell:
        "python bin/data_processing/foldx_variants.py --yaml {input.structures} {input.pdb} > {output.muts} 2> {log}"

checkpoint foldx_split:
    """
    Split variants lists into subsections to parralelise FoldX
    """
    input:
        "data/foldx/{gene}/individual_list"

    output:
        directory("data/foldx/{gene}/processing")

    params:
        n_lines = config['foldx']['variants_per_run']

    log:
        "logs/foldx_split/{gene}.log"

    shell:
        """
        mkdir data/foldx/{wildcards.gene}/processing &> {log}
        split -l {params.n_lines} data/foldx/{wildcards.gene}/individual_list data/foldx/{wildcards.gene}/processing/individual_list_ &> {log}
        """

rule foldx_repair:
    """
    Run FoldX RepairPDB command on PDB files
    """
    input:
        pdb="data/pdb/{gene}.pdb"

    output:
        "data/foldx/{gene}/{gene}_Repair.pdb",
        "data/foldx/{gene}/{gene}_Repair.fxout"

    resources:
        mem_mb = 8000

    log:
        "logs/foldx_repair/{gene}.log"

    shell:
        "foldx --command=RepairPDB --pdb={wildcards.gene}.pdb --pdb-dir=data/pdb --clean-mode=3 --output-dir=data/foldx/{wildcards.gene} &> {log}"

rule foldx_model:
    """
    Run FoldX BuildModel on a PDB and a paired list of variants
    """
    input:
        pdb="data/foldx/{gene}/{gene}_Repair.pdb",
        muts="data/foldx/{gene}/processing/individual_list_{n}"

    output:
        "data/foldx/{gene}/processing/Average_{n}_{gene}_Repair.fxout",
        "data/foldx/{gene}/processing/Dif_{n}_{gene}_Repair.fxout",
        "data/foldx/{gene}/processing/Raw_{n}_{gene}_Repair.fxout",
        "data/foldx/{gene}/processing/PdbList_{n}_{gene}_Repair.fxout"

    resources:
        mem_mb = 16000

    log:
        "logs/foldx_model/{gene}_{n}.log"

    shell:
        'foldx --command=BuildModel --pdb={wildcards.gene}_Repair.pdb --pdb-dir=data/foldx/{wildcards.gene} --mutant-file={input.muts} --output-file="{wildcards.n}" --output-dir=data/foldx/{wildcards.gene}/processing --numberOfRuns=3 --clean-mode=3 --out-pdb=false &> {log}'

def get_foldx_split_files(wildcards):
    """
    Retrieve the IDs of split FoldX jobs
    """
    checkpoint_outdir = checkpoints.foldx_split.get(gene=wildcards.gene).output[0]
    fx_output = expand('data/foldx/{gene}/processing/{fi}_{n}_{gene}_Repair.fxout',
                    gene=wildcards.gene,
                    n=glob_wildcards(os.path.join(checkpoint_outdir, "individual_list_{n}")).n,
                    fi=('Average', 'Dif', 'Raw'))
    in_lists = expand('data/foldx/{gene}/processing/individual_list_{n}',
                    gene=wildcards.gene,
                    n=glob_wildcards(os.path.join(checkpoint_outdir, "individual_list_{n}")).n)
    return fx_output + in_lists

rule foldx_combine:
    """
    Combined output of split FoldX model results for a gene
    """
    input:
        get_foldx_split_files

    output:
        "data/foldx/{gene}/average_{gene}.fxout",
        "data/foldx/{gene}/dif_{gene}.fxout",
        "data/foldx/{gene}/raw_{gene}.fxout"

    log:
        "logs/foldx_combine/{gene}.log"

    shell:
        """
        python bin/data_processing/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Average_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* > data/foldx/{wildcards.gene}/average_{wildcards.gene}.fxout 2> {log}

        python bin/data_processing/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Dif_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* --type=dif > data/foldx/{wildcards.gene}/dif_{wildcards.gene}.fxout 2> {log}

        python bin/data_processing/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Raw_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* --type=raw > data/foldx/{wildcards.gene}/raw_{wildcards.gene}.fxout 2> {log}
        """
