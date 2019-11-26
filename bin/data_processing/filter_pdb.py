#!/usr/bin/env python3
"""
Filter PDB file based on section list
"""
import sys
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from subtypes_utils import SectionSelecter, import_sections

def main(args):
    """Main script"""
    pdb_name = Path(args.pdb).stem
    # deal with FoldX repaired PDBs
    if pdb_name.endswith('_Repair'):
        pdb_name = pdb_name.replace('_Repair', '')

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_name, args.pdb)

    sections = import_sections(args.yaml, pdb_name)

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    pdbio.save(sys.stdout, select=SectionSelecter(sections))

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file")

    parser.add_argument('--yaml', '-y',
                        help=("YAML file detailing regions to process for a set of genes or "
                              "these sections in raw YAML strings"))

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
