#!/usr/bin/env python3
"""
Project various mutational landscape metrics onto protein structures
"""
import os
import argparse
from pathlib import Path
from io import StringIO

import pandas as pd
import numpy as np
import pymol2
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from subtypes_utils import SectionSelecter, import_sections, gene_to_filename
from colour_spectrum import ColourSpectrum

def main(args):
    """Main script"""
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass # Don't mind if it already exists

    dms = pd.read_csv(args.data, sep='\t')
    colourer = ColourSpectrum.linear(min(dms[args.property]), max(dms[args.property]),
                                     midpoint=0, colours='RdBu', na_colour='0xC0C0C0')

    if args.scale:
        pass

    else:
        pdb_name = args.gene or Path(args.pdb).stem

        # deal with FoldX repaired PDBs
        if pdb_name.endswith('_Repair'):
            pdb_name = pdb_name.replace('_Repair', '')

        sections = import_sections(args.structure_yaml, pdb_name)
        pdb_string = get_pdb_string(args.pdb, select=SectionSelecter(sections))

        dms = dms[[gene_to_filename(x) == pdb_name for x in dms.gene]]
        dms['pdb_position'], dms['pdb_chain'] = position_offsetter(dms.position, sections)
        dms = dms.dropna(subset=['pdb_position'])

        project_spectrum(f'{args.output_dir}/{pdb_name}_pc1.png', pdb_string, dms.pdb_chain,
                         dms.pdb_position, dms.PC1, colourer)

def get_pdb_string(pdb_path, select=None):
    """
    Get string representation of a PDB file, filtered using select as in Bio.PDB.PDBIO
    """
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('_', pdb_path)

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    with StringIO('PDB') as virtual_pdb_file:
        pdbio.save(virtual_pdb_file, select=select)
        return virtual_pdb_file.getvalue()

def position_offsetter(uniprot_positions, sections):
    """
    Convert Uniprot positions to PDB positions given a list of sections determining the
    offsets. Each section dict must have a chain, integer offset (PDB -> uniprot) and
    two item region item.
    """
    pdb_position = []
    chain = []
    for pos in uniprot_positions:
        for sec in sections:
            reg = sec['region'] if 'region' in sec else [float('-inf'), float('inf')]
            if reg[0] <= pos - sec['offset'] <= reg[1]:
                pdb_position.append(pos - sec['offset'])
                chain.append(sec['chain'])
                break
        else:
            pdb_position.append(np.nan)
            chain.append(np.nan)
    return pdb_position, chain

def project_spectrum(output_file, pdb_string, chains, positions, values, colourer):
    """
    Project a given colouring onto a PDB file via PyMOL
    """
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(pdb_string, 'prot')
        # give a values outside the range to get the NA colour - a bit hacky
        pymol.cmd.color(colourer(max(values) + 1), 'prot')

        for chn, pos, val in zip(chains, positions, values):
            pymol.cmd.color(colourer(val), f'prot and chain {chn} and resi {int(pos)}')

        pymol.cmd.png(output_file, 1000, 800, dpi=150, ray=1)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file")

    parser.add_argument('property', metavar='R', help="Property to project onto the structure")

    parser.add_argument('--scale', '-s', action='store_true',
                        help=("Output the universal colour scale for this property "
                              "instead of a projection"))

    parser.add_argument('--gene', '-g', help="Gene name (use the PDB file prefix by default)")

    parser.add_argument('--structure_yaml', '-y',
                        help="YAML config file describing structure regions",
                        default='meta/structures.yaml')

    parser.add_argument('--data', '-d', help="Landscape data",
                        default='data/combined_mutational_scans.tsv')

    parser.add_argument('--output_dir', '-o', help="Directory to output PNG files",
                        default='.')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
