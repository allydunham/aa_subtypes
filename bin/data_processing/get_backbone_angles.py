#!/usr/bin/env python3
"""
Extract backbone angles from PDB files
The config YAML should have dicts matching filenames, containing a sections list
giving the chain, offset and region of each section of the PDB file to process.
These details are the equivalent to given in the --raw argument, where the sections
list should be given in YAML format: "[{chain: A, offset: 0, region: [0, 10]}, ...]",
where region is optional, corresponds to PDB numbering and defines a closed interval.
"""
import argparse
import sys
from math import pi
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from ruamel.yaml import YAML

RAD_FACTOR = 180/pi
HEADER = 'chain\tpdb_position\taa\tposition\tphi\tpsi'

def write_backbone_angles(chain, region=None, offset=0,
                          outfile=sys.stdout, header=False):
    """
    Write Psi/Phi angles from a pdb file
    """
    if region is None:
        region = (0, float('inf'))

    polypeptide_builder = PPBuilder()
    polypeptides = polypeptide_builder.build_peptides(chain)

    if header:
        print(HEADER, file=outfile)

    for peptide in polypeptides:
        angles = peptide.get_phi_psi_list()
        for residue, (phi, psi) in zip(peptide, angles):
            position = residue.get_id()[1]
            if position >= region[0] and position <= region[1]:
                print(chain.id, position,
                      seq1(residue.get_resname()),
                      position + offset,
                      'NA' if phi is None else phi * RAD_FACTOR,
                      'NA' if psi is None else psi * RAD_FACTOR,
                      sep='\t', file=outfile)

def main(args):
    """Main script"""
    yaml = YAML(typ='safe')
    pdb_parser = PDBParser()

    pdb_name = Path(args.pdb).stem
    # deal with FoldX repaired PDBs
    if pdb_name.endswith('_Repair'):
        pdb_name = pdb_name.replace('_Repair', '')

    structure = pdb_parser.get_structure(pdb_name, args.pdb)

    sections = None
    if args.yaml:
        with open(args.yaml, 'r') as yaml_file:
            sections = yaml.load(yaml_file)[pdb_name]['sections']
    elif args.raw:
        sections = yaml.load(args.raw)

    if sections is None:
        print(f'# Backbone angles for {args.pdb}', file=sys.stdout)
        print(HEADER, file=sys.stdout)
        for chain in structure[0]:
            write_backbone_angles(chain)
    else:
        print(f'# Backbone angles for {args.pdb} in the following regions:')
        for section in sections:
            region = section['region'] if 'region' in section else 'None'
            print(f'# Chain={section["chain"]}, Offset={section["offset"]}, Region={region}',
                  file=sys.stdout)

        print(HEADER, file=sys.stdout)
        for section in sections:
            region = section['region'] if 'region' in section else None
            write_backbone_angles(structure[0][section['chain']],
                                  region=region, offset=section['offset'])


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file(s)")

    parser.add_argument('--yaml', '-y',
                        help="YAML file detailing regions to process for a set of genes")

    parser.add_argument('--raw', '-r',
                        help="Raw YAML input detailing the regions for the input PDB")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
