#!/usr/bin/env python3
"""
Filter PDB file based on section list
"""
import sys
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import Select, PDBIO
from ruamel.yaml import YAML

class SectionSelecter(Select):
    """
    Selecter based on section input in the sytle of meta/structures.yaml
    """
    def __init__(self, sections):
        self.sections = sections

    def accept_residue(self, residue):
        if self.sections is None:
            return 1

        chain = residue.full_id[2]
        position = residue.id[1]
        for section in self.sections:
            if (chain == section['chain'] and
                    section['region'][0] <= position <= section['region'][1]):
                return 1
        return 0

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

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    pdbio.save(sys.stdout, select=SectionSelecter(sections))

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file")

    parser.add_argument('--yaml', '-y',
                        help="YAML file detailing regions to process for a set of genes")

    parser.add_argument('--raw', '-r',
                        help="Raw YAML input detailing the regions for the input PDB")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
