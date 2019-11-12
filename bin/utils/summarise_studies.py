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
    gene_data = {}
    for study_dir in args.dirs:
        study_dir = study_dir.rstrip('/')
        study = study_dir.split('/')[-1]

        with open(f'{study_dir}/{study}.yaml', 'r') as yaml_file:
            meta = yaml.load(yaml_file)

        gene = meta['gene']

        # calculate study data
        variants = pd.read_csv(f'{study_dir}/{study}.tsv', sep='\t')
        positions = variants.position.nunique()
        completeness = round(positions/len(meta['seq']), 3)
        total = variants[variants['class'] == 'Missense'].shape[0]
        coverage = round(total/(positions*19), 3)
        nonsense = 'Nonsense' in variants['class'].values
        synonymous = 'Synonymous' in variants['class'].values

        study_data['Study'].append(meta['study'])
        study_data['Authour'].append(meta['authour'])
        study_data['Year'].append(meta['year'])
        study_data['Gene'].append(gene)
        study_data['Function'].append(meta['gene_type'])
        study_data['Species'].append(meta['species'])
        study_data['Experiment'].append(meta['experiment'])
        study_data['Transform'].append(meta['transform'])
        study_data['Filtered'].append(meta['qc']['filter'])
        study_data['Gene Length'].append(len(meta['seq']))
        study_data['Mutated Positions'].append(positions)
        study_data['Completeness (Missense)'].append(completeness)
        study_data['Variants (Missense)'].append(total)
        study_data['Coverage (Missense)'].append(coverage)
        study_data['Nonsense'].append(nonsense)
        study_data['Synonymous'].append(synonymous)

        # Add to gene data
        if not meta['qc']['filter'] and gene in gene_data:
            gene_data[gene]['df'] = pd.merge(gene_data[gene]['df'],
                                             variants[['position', 'wt', 'mut', 'class']],
                                             how='outer', on=['position', 'wt', 'mut', 'class'])

            gene_data[gene]['studies'] += 1

            if not gene_data[gene]['length'] == len(meta['seq']):
                raise ValueError(f'{gene}: Lengths do not match')

            if not gene_data[gene]['function'] == meta['gene_type']:
                raise ValueError(f'{gene}: Functions do not match')

            if not gene_data[gene]['species'] == meta['species']:
                raise ValueError(f'{gene}: Species do not match')

        elif not meta['qc']['filter']:
            gene_data[gene] = {'df': variants[['position', 'wt', 'mut', 'class']],
                               'studies': 1, 'length': len(meta['seq']),
                               'uniprot_id': meta['uniprot_id'],
                               'function': meta['gene_type'],
                               'species': meta['species']}

    # Write study table
    study_data = pd.DataFrame(study_data)
    study_data.to_csv(args.studies, sep='\t', index=False)

    # Calculate gene table
    gene_df = defaultdict(list)
    for gene, data in gene_data.items():
        positions = data['df'].position.nunique()
        variants = data['df'][~(data['df'].mut == '*')].shape[0]
        completeness = round(positions / data['length'], 3)
        coverage = round(data['df'][data['df']['class'] == 'Missense'].shape[0]/(positions*19), 3)

        gene_df['Gene'].append(gene)
        gene_df['Uniprot ID'].append(data['uniprot_id'])
        gene_df['Function'].append(data['function'])
        gene_df['Species'].append(data['species'])
        gene_df['Studies'].append(data['studies'])
        gene_df['Length'].append(data['length'])
        gene_df['Mutated Positions'].append(positions)
        gene_df['Completeness (Missense)'].append(completeness)
        gene_df['Variants (Missense)'].append(variants)
        gene_df['Coverage (Missense)'].append(coverage)

    gene_df = pd.DataFrame(gene_df)
    gene_df.to_csv(args.genes, sep='\t', index=False)

    # Calculate Summary
    with open(args.summary, 'w') as summary_file:
        print('Summary statistics for Subtypes project', file=summary_file)
        print('Total Studies:', study_data.shape[0], file=summary_file)
        print('Unfiltered Studies:', sum(~study_data['Filtered']), file=summary_file)
        print('Species:', gene_df.Species.nunique(),
              f"({', '.join(set(gene_df.Species))})", file=summary_file)
        print('Genes:', gene_df.shape[0], file=summary_file)
        print('Positions:', sum(gene_df['Mutated Positions']), file=summary_file)
        print('Missense Variants:', sum(gene_df['Variants (Missense)']), file=summary_file)
        print('Average Completeness:',
              sum(gene_df['Completeness (Missense)'])/gene_df.shape[0],
              file=summary_file)
        print('Average Coverage:',
              sum(gene_df['Coverage (Missense)'])/gene_df.shape[0],
              file=summary_file)


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dirs', metavar='D', nargs='+', help="Input study directories")

    parser.add_argument('--studies', '-s', help="Output file for study table",
                        default='studies.tsv')

    parser.add_argument('--genes', '-g', help="Output file for genes table",
                        default='genes.tsv')

    parser.add_argument('--summary', '-u', help="Output file for summary table",
                        default='summary.tsv')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
