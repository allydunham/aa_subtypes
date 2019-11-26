#!/usr/bin/env python3
"""
Template Script
"""
import sys
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from subtypes_utils import import_sections

AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'

def main(args):
    """Main script"""
    pdb_parser = PDBParser()

    pdb_name = Path(args.pdb).stem
    # deal with FoldX repaired PDBs
    if pdb_name.endswith('_Repair'):
        pdb_name = pdb_name.replace('_Repair', '')

    structure = pdb_parser.get_structure(pdb_name, args.pdb)

    sections = import_sections(args.yaml, pdb_name)

    variants = []
    if sections is not None:
        for section in sections:
            filter_region = 'region' in section
            for residue in structure[0][section['chain']]:
                if not residue.id[0] == ' ':
                    continue # Filter HETATMs

                position = int(residue.id[1])
                amino_acid = seq1(residue.get_resname())

                if (filter_region and
                        (position > section['region'][1] or
                         position < section['region'][0])):
                    continue

                variants.extend([f"{amino_acid}{section['chain']}{position}{x}" for
                                 x in AA_ALPHABET if not x == amino_acid])
    else:
        for chain in structure[0]:
            for residue in chain:
                if not residue.id[0] == ' ':
                    continue # Filter HETATMs

                position = int(residue.id[1])
                amino_acid = seq1(residue.get_resname())

                variants.extend([f"{amino_acid}{chain.id}{position}{x}" for
                                 x in AA_ALPHABET if not x == amino_acid])

    print(*variants, sep=';\n', end=';\n', file=sys.stdout)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB")

    parser.add_argument('--yaml', '-y', help="YAML file/raw YAML defining PDB sections to consider")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
