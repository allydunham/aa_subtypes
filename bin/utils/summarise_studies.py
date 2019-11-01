#!/usr/bin/env python3
"""
Script generating summary tables of the different studies used in the project
"""
import argparse
from collections import defaultdict

import pandas as pd
from ruamel.yaml import YAML

def main(args):
    """Main script"""
    yaml = YAML(typ='safe')

    study_data = defaultdict(list)
    for study_dir in args.dirs:
        study_dir = study_dir.rstrip('/')
        study = study_dir.split('/')[-1]

        with open(f'{study_dir}/{study}.yaml', 'r') as yaml_file:
            meta = yaml.load(yaml_file)

        variants = pd.read_csv(f'{study_dir}/{study}.tsv', sep='\t')
        positions = 0
        coverage = 0
        completeness = 0

        study_data['Study'].append(meta['study'])
        study_data['Authour'].append(meta['authour'])
        study_data['Year'].append(meta['year'])
        study_data['Gene'].append(meta['gene'])
        study_data['Function'].append(meta['gene_type'])
        study_data['Species'].append(meta['species'])
        study_data['Experiment'].append(meta['experiment'])
        study_data['Transform'].append(meta['transform'])
        study_data['Filtered'].append(meta['qc']['filter'])
        study_data['Positions'].append(positions)
        study_data['Coverage'].append(coverage)
        study_data['Completeness'].append(completeness)

    study_data = pd.DataFrame(study_data)
    study_data.to_csv(args.studies, sep='\t', index=False)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dirs', metavar='D', nargs='+', help="Input study directories")

    parser.add_argument('--studies', '-s', help="Output file for study table",
                        default='studies.tsv')

    parser.add_argument('--genes', '-g', help="Output file for genes table [NOT IMPLEMENTED]",
                        default='genes.tsv')

    parser.add_argument('--summary', '-u', help="Output file for summary table",
                        default='summary.tsv')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
