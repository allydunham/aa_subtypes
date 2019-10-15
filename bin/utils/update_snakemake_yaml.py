#!/usr/bin/env python3
"""
Script to update the snakemake config YAML with all the studies/genes/structures present
"""
import argparse
import os
from ruamel.yaml import YAML

def main(args):
    """Main script"""
    yaml = YAML(typ='safe')

    try:
        with open(args.yaml, 'r') as yaml_file:
            conf = yaml.load(yaml_file)

        if conf is None:
            conf = {}

    except FileNotFoundError:
        conf = {}

    conf['studies'] = []
    conf['genes'] = []

    for study in os.listdir(args.studies):
        with open(f'{args.studies}/{study}/{study}.yaml', 'r') as yaml_file:
            meta = yaml.load(yaml_file)

        conf['studies'].append(meta['study'])
        conf['genes'].append(meta['gene'].lower())

    with open(args.yaml, 'w') as yaml_file:
        yaml.dump(conf, yaml_file)


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-y', '--yaml', default='snakemake.yaml',
                        help="Current snakemake config YAML")

    parser.add_argument('-s', '--studies', default='data/studies',
                        help='Root directory containing study folders')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
