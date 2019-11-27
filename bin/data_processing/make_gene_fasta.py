#!/usr/bin/env python3
"""
Generate a fasta file for a gene, given a set of study yamls
"""
import argparse
from ruamel.yaml import YAML

from subtypes_utils import gene_to_filename

FASTA_LINE_LENGTH = 80

def main(args):
    """Main script"""
    yaml = YAML(typ='safe')
    seq = None
    gene = None
    for study_yaml in args.yaml:
        with open(study_yaml, 'r') as yaml_file:
            conf = yaml.load(yaml_file)

        if seq is None:
            seq = conf['seq']
            gene = conf['gene']
        elif not gene == conf['gene']:
            raise ValueError(f"Studies are for different genes")
        elif not seq == conf['seq']:
            raise ValueError(f"Studies have different sequences for {gene}")

    print(f">{gene_to_filename(gene)}")
    for i in range(0, len(seq), FASTA_LINE_LENGTH):
        print(seq[i:(i + FASTA_LINE_LENGTH)])

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('yaml', metavar='Y', nargs='+', help="Input study config YAML file(s)")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
