#!/usr/bin/env python3
"""
Generate a fasta file for a DMS study from its meta YAML file
"""
import sys
import argparse
from ruamel.yaml import YAML

def main(args):
    """Main script"""
    yaml = YAML(typ='safe')

    with open(args.yaml, 'r') as yaml_file:
        meta = yaml.load(yaml_file)

    print(f">{meta['gene']}", file=sys.stdout)
    for i in range(0, len(meta['seq']), args.length):
        print(meta['seq'][i:(i + args.length)], file=sys.stdout)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('yaml', metavar='Y', help="YAML file with meta details about a study")

    parser.add_argument('-l', '--length', type=int, default=80,
                        help="Maximum number of characters per line")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
