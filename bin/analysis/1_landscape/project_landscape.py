#!/usr/bin/env python3
"""
Project mutational landscape metrics onto protein structures and/or produce the corresponding
colour scales.
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

# Give the scale and midpoint to use for various properties, default is (bwr, 0)
PROPERTIES = {
    'PC1': ('RdBu', 0), 'total_energy': ('PiYG', 0), 'mean_sift': ('PuRd', None)
}

def main(args):
    """Main script"""
    # Process arguments
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass # Don't mind if it already exists

    if not (args.colourbar or args.pdb):
        raise ValueError('Provide input PDB files or plot colourbars (--colourbar)')

    if args.gene:
        if len(args.gene) != len(args.pdb):
            raise ValueError(('--genes/-g and --pdb/-p lists must be the same length when '
                              'using --genes'))
        genes = args.genes
    elif args.pdb:
        # Cover the common case of FoldX adding _Repair
        genes = [Path(i).stem.rsplit('_Repair', 1)[0] for i in args.pdb]

    # import DMS data
    dms = pd.read_csv(args.data, sep='\t')

    # Calculate (and optionally plot) global scales
    if not args.local_scale:
        colourers = [get_colourer(l, dms) for l in args.properties]
        if args.colourbar:
            for colour_scale in colourers:
                fig, _ = colour_scale.plot()
                fig.savefig(f'{args.output_dir}/{colour_scale.name}_colourbar{args.suffix}.pdf',
                            bbox_inches='tight')

    # Loop over PDBs/Properties and project, including determining local scales if needed
    if args.pdb:
        for gene, pdb_path in zip(genes, args.pdb):
            pdb_string, pdb_dms = get_pdb_data(pdb_path, gene, args.structure_yaml, dms)
            if args.local_scale:
                colourers = {l: get_colourer(l, dms) for l in args.properties}
                if args.colourbar:
                    for colour_scale in colourers:
                        fig, _ = colour_scale.plot()
                        path = (f'{args.output_dir}/{gene}_{colour_scale.name}'
                                f'_colourbar{args.suffix}.pdf')
                        fig.savefig(path, bbox_inches='tight')

            for landscape_property in args.properties:
                colourer = colourers[landscape_property]
                project_spectrum(f'{args.output_dir}/{gene}_{landscape_property}{args.suffix}.png',
                                    pdb_string, pdb_dms.pdb_chain, pdb_dms.pdb_position,
                                    pdb_dms[landscape_property], colourer)

def get_colourer(landscape_property, dms):
    """
    Generate a ColourSpectrum.linear for a given column in the dms data, using
    the colour scheme set in PROPERTIES with the default being RdBu with midpoint=0
    """
    minimum = min(dms[landscape_property])
    maximum = max(dms[landscape_property])
    spectrum, midpoint = PROPERTIES.get(landscape_property, ('bwr', 0))
    return ColourSpectrum(minimum, maximum, midpoint=midpoint, name=landscape_property,
                          colourmap=spectrum, na_colour='0xC0C0C0')

def get_pdb_data(path, gene, section_yaml, dms_df):
    """
    Import PDB data from file and a matching deep scanning dataframe, filtering
    to the desired regions.
    """
    sections = import_sections(section_yaml, gene)
    pdb_string = get_pdb_string(path, select=SectionSelecter(sections))

    # Filter dms to only include the pdb region
    pdb_dms = dms_df[[gene_to_filename(x) == gene for x in dms_df.gene]].copy()
    pdb_dms['pdb_position'], pdb_dms['pdb_chain'] = position_offsetter(pdb_dms.position,
                                                                       sections)
    pdb_dms = pdb_dms.dropna(subset=['pdb_position'])

    return pdb_string, pdb_dms

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
        pymol.cmd.color(colourer.na_colour, 'prot')

        for chn, pos, val in zip(chains, positions, values):
            pymol.cmd.color(colourer(val), f'prot and chain {chn} and resi {int(pos)}')

        pymol.cmd.png(output_file, 1000, 800, dpi=150, ray=1)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('properties', metavar='P', nargs='+',
                         help="Property(s) to project")

    parser.add_argument('--pdb', '-p', nargs='+', help="Input PDB file(s)")

    parser.add_argument('--colourbar', '-c', action='store_true',
                        help=("Output the universal colour scale for this property."
                              "Produce one spectrum per property per PDB if using --local_scale, "
                              "otherwise one per property"))

    parser.add_argument('--local_scale', '-s', action='store_true',
                        help="Define scales on values of the property in each structure")

    parser.add_argument('--gene', '-g', nargs='+',
                        help=("Manual gene name(s) for each PDB file. If not provided"
                              "the file stem, sans FoldXs _Repair suffix is used."))

    parser.add_argument('--structure_yaml', '-y',
                        help="YAML config file describing structure regions",
                        default='meta/structures.yaml')

    parser.add_argument('--data', '-d', help="Landscape data",
                        default='data/combined_mutational_scans.tsv')

    parser.add_argument('--output_dir', '-o', help="Directory to output PNG files",
                        default='.')

    parser.add_argument('--suffix', '-u', help="Suffix to add to output plots",
                        default='')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
