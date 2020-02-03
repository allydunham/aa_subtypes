#!/usr/bin/env python3
"""
Project mutational landscape metrics onto protein structures and/or produce the corresponding
colour scales.
"""
import os
import argparse
from pathlib import Path

import pandas as pd
import pymol2
from ruamel.yaml import YAML

import subtypes_utils as su
from pymol_utils import project_landscape

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
        genes = args.gene
    elif args.pdb:
        # Cover the common case of FoldX adding _Repair
        genes = [Path(i).stem.rsplit('_Repair', 1)[0] for i in args.pdb]

    # import data
    yaml = YAML(typ='safe')
    with open(args.structure_yaml, 'r') as yaml_file:
        sections = {k: v['sections'] for k, v in yaml.load(yaml_file).items()}

    dms = pd.read_csv(args.data, sep='\t')
    dms['pdb_chain'], dms['pdb_position'] = su.dms_pdb_positions(dms, sections)
    dms = dms.dropna(subset=['pdb_position'])

    # Calculate (and optionally plot) global scales
    if not args.local_scale:
        colourers = make_colourers(dms, args.properties, suffix=args.suffix,
                                   save=args.output_dir if args.colourbar else False)

    # Loop over PDBs/Properties and project, including determining local scales if needed
    if args.pdb:
        for gene, pdb_path in zip(genes, args.pdb):
            pdb_string = su.get_pdb_string(pdb_path, select=su.SectionSelecter(sections[gene]))
            pdb_dms = dms[dms.gene.apply(su.gene_to_filename) == gene].copy()

            if args.local_scale:
                colourers = make_colourers(pdb_dms, args.properties, save=args.output_dir,
                                           suffix=args.suffix)

            for landscape_property in args.properties:
                colourer = colourers[landscape_property]
                with pymol2.PyMOL() as pymol:
                    pymol.cmd.read_pdbstr(pdb_string, 'prot')
                    project_landscape(pymol.cmd, pdb_dms.pdb_chain, pdb_dms.pdb_position,
                                      pdb_dms[landscape_property], colourer=colourer)
                    pymol.cmd.png(f'{args.output_dir}/{gene}_{landscape_property}{args.suffix}.png',
                                  1000, 800, dpi=150, ray=1)

def make_colourers(dms, properties, save=None, suffix=None):
    """
    Generate a dictionary of colourers and optionally save images of them. If saving save should
    be set to the directory to output to
    """
    colourers = {l: su.SubtypesColourSpectrum(l, dms[l]) for l in properties}
    if save:
        for name, colour_scale in colourers.items():
            fig, _ = colour_scale.plot()
            fig.savefig(f'{save}/{name}_colourbar{suffix}.pdf',
                        bbox_inches='tight')

    return colourers

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('properties', metavar='P', nargs='+',
                        help="Property(s) to project")

    parser.add_argument('--pdb', '-p', nargs='+', help="Input PDB file(s)")

    parser.add_argument('--colourbar', '-c', action='store_true',
                        help=("Output global colour scales. Has no effect if using --local_scale"
                              "where colour scale plots are always produced."))

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
