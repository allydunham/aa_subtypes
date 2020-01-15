#!/usr/bin/env python3
"""
Project various mutational landscape metrics onto protein structures
"""
import os
import argparse
from pathlib import Path
from io import StringIO
from bisect import bisect

import pandas as pd
import numpy as np
import pymol2
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from subtypes_utils import SectionSelecter, import_sections, gene_to_filename

def main(args):
    """Main script"""
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass # Don't mind if it already exists


    pdb_name = args.gene or Path(args.pdb).stem

    # deal with FoldX repaired PDBs
    if pdb_name.endswith('_Repair'):
        pdb_name = pdb_name.replace('_Repair', '')

    sections = import_sections(args.structure_config, pdb_name)

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_name, args.pdb)

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    with StringIO('PDB') as virtual_pdb_file:
        pdbio.save(virtual_pdb_file, select=SectionSelecter(sections))
        pdb_string = virtual_pdb_file.getvalue()

    dms = pd.read_csv(args.data, sep='\t')
    dms = dms[[gene_to_filename(x) == pdb_name for x in dms.gene]]
    dms['pdb_position'], dms['pdb_chain'] = position_offsetter(dms.position, sections)
    dms = dms.dropna(subset=['pdb_position'])

    colourer = get_colour_spectrum(dms.PC1, [(239, 138, 98), (247, 247, 247), (103, 169, 207)])
    project_spectrum(f'{args.output_dir}/{pdb_name}_pc1.png', pdb_string, dms.pdb_chain,
                     dms.pdb_position, dms.PC1, colourer)

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

def get_colour_spectrum(values, colours, sym=None, na_colour='0xC0C0C0'):
    """
    Generate a sensible colour_spectrum function for a given input vector.
    Sym makes the scale symetric about a value, extending as far as the most
    distant value from that point.
    """
    min_val = min(values)
    max_val = max(values)
    n_colours = len(colours)

    if sym is not None:
        if n_colours % 2 == 0:
            mid_colour = rgb_interpolate(colours[n_colours // 2 - 1],
                                         colours[n_colours // 2],
                                         0.5)
            colours.insert(n_colours // 2, mid_colour)

        diff = max(abs(max_val - sym), abs(min_val - sym))
        max_val = sym + diff
        min_val = sym - diff

    colour_div = (max_val - min_val)/n_colours
    colour_values = [min_val + colour_div * x for x in range(n_colours)]
    return colour_spectrum(colour_values, colours, outlier_colour=na_colour)


def colour_spectrum(values, colours, outlier_colour='0xC0C0C0'):
    """
    Create a function converting numeric values to hex RGB codes.

    values: iterable of numeric values
    colours: iterable of RGB tuples corresponding to values
    outlier_colour: RGB tuple or Hex code of colour to return for out out range values
    """
    if len(values) != len(colours):
        raise ValueError('values and colours must be the same length')

    if not isinstance(outlier_colour, str):
        outlier_colour = rgb_to_hex(outlier_colour)

    values = sorted(list(zip(values, colours)), key=lambda x: x[0])
    values, colours = zip(*values)

    def spectrum(val):
        if val in values:
            return rgb_to_hex(colours[values.index(val)])

        ind = bisect(values, val)

        if ind == 0 or ind >= len(values):
            return outlier_colour

        prop = (val - values[ind - 1]) / (values[ind] - values[ind - 1])
        res = rgb_interpolate(colours[ind - 1], colours[ind], prop)

        return rgb_to_hex(res)

    return spectrum

def rgb_interpolate(low, high, prop):
    """
    Interpolate between two RGB tuples
    """
    return [rgb_clamp(y + (x - y) * prop) for x, y in zip(low, high)]


def rgb_clamp(colour_value):
    """
    Clamp a value to integers on the RGB 0-255 range
    """
    return int(min(255, max(0, colour_value)))

def rgb_to_hex(rgb):
    """
    Covert RGB tuple to Hex code
    """
    return f'0x{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}'

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Input PDB file")

    parser.add_argument('--gene', '-g', help="Gene name (use the PDB file prefix by default)")

    parser.add_argument('--structure_config', '-s',
                        help="YAML config file describing structure regions",
                        default='meta/structures.yaml')

    parser.add_argument('--data', '-d', help="Landscape data",
                        default='data/combined_mutational_scans.tsv')

    parser.add_argument('--output_dir', '-o', help="Directory to output PNG files",
                        default='figures/1_landscape/pdb')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
