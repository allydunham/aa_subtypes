#!/usr/bin/env python3
"""
Calculate chemical environment profiles for all positions in a PDB file.
Can select sections to use using the format of meta/structures.yaml, although
offsets are not applied so as to keep numbering consistent with PDBs so
this must be applied separately to map to Uniprot sequences.
"""
import sys
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.Alphabet.IUPAC import protein as protein_alphabet
from ruamel.yaml import YAML

import chemical_environment as ce
from filter_pdb import SectionSelecter

# TODO loading bars/progress stuff

def main(args):
    """Main script"""
    pdb_name = Path(args.pdb).stem
    # deal with FoldX repaired PDBs
    if pdb_name.endswith('_Repair'):
        pdb_name = pdb_name.replace('_Repair', '')

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_name, args.pdb)

    sections = None
    yaml = YAML(typ='safe')
    if args.yaml:
        with open(args.yaml, 'r') as yaml_file:
            sections = yaml.load(yaml_file)[pdb_name]['sections']
    elif args.raw:
        sections = yaml.load(args.raw)

    selecter = SectionSelecter(sections, drop_hetero=True)
    residues = [r for c in structure[0] for r in c if selecter.accept_residue(r)]

    if args.k_nearest:
        prof_func = lambda x: ce.k_nearest_residues(x, k=args.k_nearest)
        profile_cols = [f'nearest_{args.k_nearest}_{aa}' for aa in protein_alphabet.letters]
    elif args.angstroms:
        prof_func = lambda x: ce.within_distance(x, max_dist=args.angstroms)
        profile_cols = [f'within_{args.angstroms}_{aa}'.replace('.', '_') for aa in
                        protein_alphabet.letters]

    print('chain', 'position', 'wt', *profile_cols, sep='\t', file=sys.stdout)
    for residue, profile in zip(residues, prof_func(residues)):
        print(residue.full_id[2], residue.id[1], seq1(residue.get_resname()),
              *profile, sep='\t', file=sys.stdout)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file")

    parser.add_argument('--yaml', '-y',
                        help="YAML file detailing regions to process for a set of genes")

    parser.add_argument('--raw', '-r',
                        help="Raw YAML input detailing the regions for the input PDB")

    parser.add_argument('--k_nearest', '-k', type=int,
                        help='Profile based on k nearest residues')

    parser.add_argument('--angstroms', '-a', type=float,
                        help='Profile based on count of residues within a angstroms')

    args = parser.parse_args()

    if not (args.k_nearest or args.angstroms):
        raise ValueError('Please select a profile method (--k_nearest or --angstroms)')

    elif args.k_nearest and args.angstroms:
        raise ValueError('Use either --k_nearest or --angstroms, not both')

    return args

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
