#!/usr/bin/env python3
"""
Template Script
"""
import sys
import argparse
import pandas as pd
from numpy import repeat, tile, arange

def main(args):
    """Main script"""
    if not len(args.foldx) == len(args.variants):
        raise ValueError('Each input table must be paired with an individual_list of variants')

    data_frames = []
    for df_path, variant_path in zip(args.foldx, args.variants):
        frame = pd.read_csv(df_path, sep='\t', skiprows=8)
        frame = frame.rename(lambda x: x.lower().replace(' ', '_'), axis='columns')

        with open(variant_path, 'r') as variant_file:
            variants = [x.rstrip(';\n') for x in variant_file.readlines()]

        if args.type == 'average':
            frame = frame.drop('pdb', axis='columns')
            frame['variant'] = variants
            col_order = ['chain', 'position', 'wt', 'mut'] + frame.columns.tolist()[:-1]

        elif args.type == 'dif':
            n_reps = get_fx_reps(frame.pdb.values)
            frame = frame.drop('pdb', axis='columns')
            frame['variant'] = repeat(variants, n_reps)
            col_order = ['chain', 'position', 'wt', 'mut'] + frame.columns.tolist()[:-1]

        elif args.type == 'raw':
            n_reps = get_fx_reps(frame.pdb.values)
            frame = frame.drop('pdb', axis='columns')
            frame['variant'] = repeat(variants, 2*n_reps)
            frame['source'] = tile(['mut']*n_reps + ['wt']*n_reps, frame.shape[0]//(2*n_reps))
            frame['rep'] = tile(arange(0, n_reps), frame.shape[0]//n_reps)
            col_order = ['source', 'rep', 'chain', 'position', 'wt', 'mut'] + frame.columns.tolist()[:-3]

        regex = r"(?P<wt>[A-Z])(?P<chain>.)(?P<position>[0-9]+)(?P<mut>[A-Z])"
        frame = pd.concat((frame, frame.variant.str.extract(regex)), axis=1)
        frame = frame.drop('variant', axis='columns')
        frame = frame[col_order]
        data_frames.append(frame)

    data_frame = pd.concat(data_frames)
    data_frame.to_csv(sys.stdout, sep='\t', index=False, float_format='%.6g')

def get_fx_reps(pdb):
    """
    Determine the number of FoldX replicates run from PDB column output
    """
    count = 0
    for index in pdb:
        index = index.split('_')
        if index[0] == 'WT' or not index[-2] == '1':
            # WT catches end of first PDB in Raw files
            # x[-2] catches end in Dif/Average files
            break
        count += 1
    return count

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--foldx', '-f', nargs='+', help="Files to combine", required=True)
    parser.add_argument('--variants', '-v', nargs='+', help="Variant lists for each input file",
                        required=True)

    parser.add_argument('--type', '-t', choices=['average', 'raw', 'dif'],
                        help="Type of FoldX file to process", default='average')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
