#!/usr/bin/env python3
"""
Calculate chemical environment profiles for all positions in a PDB file(s) or PDB files
linked from dm files. If a single PDB file is given output will be to stdout unless a file
suffix is supplied. Otherwise output files are written in the directories of input
PDB files.
"""
import sys
import os
import argparse
import Bio.PDB
from Bio.SeqUtils import seq1
import chemical_environment as ce
import deep_mut_tools as dmt

# TODO don't assume model 0 always?
# TODO Add keep_hetero option?
# TODO TEST!
# TODO loading bars/progress stuff

def main(args):
    """Main script"""
    if not (args.k_nearest or args.angstroms):
        raise ValueError(('No profile method selected, at least one of'
                          '--k_nearest/--angstroms must have at least one value'))

    # construct list of PDB files to process
    if args.dm:
        pdb_files = []
        for dm_file in args.input:
            root = '/'.join(dm_file.split('/')[:-1])
            dm_header = dmt.read_deep_mut_header(dm_file)

            if dm_header['pdb_id'] is None:
                print(f'Warning: No PDB IDs found in dm file {dm_file}', file=sys.stderr)
                continue

            try:
                pdb_ids = [dm_header['pdb_id'].split(':')]
            except AttributeError:
                pdb_ids = [i.split(':') for i in dm_header['pdb_id']]

            for pdb in pdb_ids:
                pdb_files.append({'pdb_file': f'{root}/{pdb[0]}/{pdb[0]}.pdb',
                                  'groups': 'all' if args.combine_chains else [[pdb[1]]],
                                  'offset': {pdb[1]:int(pdb[2])} if len(pdb) > 2 else 0})

    else:
        pdb_files = [{'pdb_file': f,
                      'groups': 'all' if args.combine_chains else 'chains',
                      'offset': 0} for f in args.input]

    # Determine output method
    write_stdout = not args.suffix and len(pdb_files) == 1
    suffix = args.suffix or '.chem_env'

    # Calculate profiles and write tables for each PDB file
    pdb_parser = Bio.PDB.PDBParser()
    for pdb in pdb_files:
        pdb_id = pdb['pdb_file'].split('/')[-1].split('.')[0]
        structure = pdb_parser.get_structure(pdb_id,
                                             pdb['pdb_file'])
        chem_envs = ce.ChemicalEnvironment(structure[0], pdb_id=pdb_id)
        chem_envs.select_residues(groups=pdb['groups'], drop_hetero=True) #args.keep_hetero

        # Generate profiles
        if args.k_nearest:
            for k in args.k_nearest:
                chem_envs.calc_k_nearest_profiles(k)

        if args.angstroms:
            for max_dist in args.angstroms:
                chem_envs.calc_within_distance_profiles(max_dist)

        if write_stdout:
            chem_envs.write_profile_table(sys.stdout, offset=pdb['offset'])

        else:
            with open(f"{os.path.splitext(pdb['pdb_file'])[0]}{suffix}", 'w') as outfile:
                chem_envs.write_profile_table(outfile, offset=pdb['offset'])

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', metavar='I', nargs='+', help="Input file(s)")

    parser.add_argument('--dm', '-d', action='store_true', help="Assume input are dm files")

    parser.add_argument('--combine_chains', '-c', action='store_true',
                        help="Process all chains together from PDB input")

    # parser.add_argument('--keep_hetero', '-k', action='store_true',
    #                     help="Include hetero atoms in profiles")

    parser.add_argument('--suffix', '-s', default='',
                        help="Suffix to add to replace .pdb with in generated files")

    parser.add_argument('--k_nearest', '-k', nargs='+', type=int,
                        help='list of k values to use for k-nearest amino acids profiles')

    parser.add_argument('--angstroms', '-a', nargs='+', type=float,
                        help='List of maximum distances (in A) to use in within distance profiles')

    parser.add_argument('-X', action='store_true',
                        help='Separator to mark start of file list if needed with list arguments')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
