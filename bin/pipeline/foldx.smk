"""
Rules for FoldX pipeline

Expects two global variables:
- AA_ALPHABET: string of the 20 AA letters) and
- STRUCTURES: dictionary of meta/structures.yaml
"""

rule foldx_variants:
    input:
        pdb="data/foldx/{gene}/{gene}.pdb",
        structures=ancient("meta/structures.yaml") # Imported generally

    output:
        muts="data/foldx/{gene}/individual_list"

    run:
        parser = PDBParser()
        structure = parser.get_structure(f'{wildcards.gene}', input.pdb)

        variants = []
        for section in STRUCTURES[wildcards.gene]['sections']:
            filter_region = 'region' in section
            for residue in structure[0][section['chain']]:
                if not residue.id[0] == ' ':
                    continue # Filter HETATMs

                position = int(residue.id[1])
                aa = seq1(residue.get_resname())

                if (filter_region and
                    (position > section['region'][1] or
                     position < section['region'][0])):
                    continue

                variants.extend([f"{aa}{section['chain']}{position}{x}" for
                                 x in AA_ALPHABET if not x == aa])

        with open(output.muts, 'w') as outfile:
            print(*variants, sep=';\n', end=';\n', file=outfile)

checkpoint foldx_split:
    input:
        "data/foldx/{gene}/individual_list"

    output:
        temp(directory("data/foldx/{gene}/processing"))

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

    resources:
        mem_mb = 4000

    log:
        "logs/foldx_repair/{gene}.log"

    shell:
        "foldx --command=RepairPDB --pdb={wildcards.gene}.pdb --pdb-dir=data/foldx/{wildcards.gene} --clean-mode=3 --output-dir=data/foldx/{wildcards.gene} &> {log}"

rule foldx_model:
    input:
        pdb="data/foldx/{gene}/{gene}_Repair.pdb",
        muts="data/foldx/{gene}/processing/individual_list_{n}"

    output:
        temp("data/foldx/{gene}/processing/Average_{n}_{gene}_Repair.fxout"),
        temp("data/foldx/{gene}/processing/Dif_{n}_{gene}_Repair.fxout"),
        temp("data/foldx/{gene}/processing/Raw_{n}_{gene}_Repair.fxout"),
        temp("data/foldx/{gene}/processing/PdbList_{n}_{gene}_Repair.fxout")

    resources:
        mem_mb = 4000

    log:
        "logs/foldx_model/{gene}_{n}.log"

    shell:
        'foldx --command=BuildModel --pdb={wildcards.gene}_Repair.pdb --pdb-dir=data/foldx/{wildcards.gene} --mutant-file={input.muts} --output-file="{wildcards.n}" --output-dir=data/foldx/{wildcards.gene}/processing --numberOfRuns=3 --clean-mode=3 --out-pdb=false &> {log}'

def get_foldx_split_files(wildcards):
    """
    Retrieve the split FoldX jobs
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
    input:
        get_foldx_split_files

    output:
        "data/foldx/{gene}/average_{gene}.fxout",
        "data/foldx/{gene}/dif_{gene}.fxout",
        "data/foldx/{gene}/raw_{gene}.fxout"

    shell:
        """
        python bin/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Average_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* > data/foldx/{wildcards.gene}/average_{wildcards.gene}.fxout

        python bin/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Dif_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* --type=dif > data/foldx/{wildcards.gene}/dif_{wildcards.gene}.fxout

        python bin/foldx_combine.py --foldx data/foldx/{wildcards.gene}/processing/Raw_*_{wildcards.gene}_Repair.fxout --variants data/foldx/{wildcards.gene}/processing/individual_list_* --type=raw > data/foldx/{wildcards.gene}/raw_{wildcards.gene}.fxout
        """
