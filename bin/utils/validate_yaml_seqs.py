#!/usr/bin/env python3
"""
Check sequences in study config YAML files (assumed to be in {study}/{study}.yaml format)
against a folder of raw fasta sequences.
Fasta files are assumed to be named in the format {species e.g. h_sapiens}_{gene}.fa/fasta
except for the special
"""
import argparse
import os

import ruamel.yaml
from Bio import SeqIO
from colorama import Fore
from colorama import Style

PROJECT_ROOT = '/Users/ally/phd/subtypes'
DEFAULT_STUDIES = f'{PROJECT_ROOT}/data/studies'
DEFAULT_FASTA = f'{PROJECT_ROOT}/meta/fasta'

def import_single_seq_fasta_dir(directory):
    """Import a directory of single sequence fasta files into a dict, keyed by the file names"""
    seqs = {}
    for fasta in os.listdir(directory):
        seq = list(SeqIO.parse(f'{directory}/{fasta}', 'fasta', ))

        if not seq:
            print(f'Warning: non-fasta formatted file found in directory ({fasta})')
            continue

        if len(seq) > 1:
            print(f'Warning: More than 1 sequence found in {fasta}, only using the first')

        name = os.path.splitext(fasta)[0]
        seqs[name] = seq[0]
    return seqs

def main(args):
    """Main script"""
    seqs = import_single_seq_fasta_dir(args.fasta)

    yaml_loader = ruamel.yaml.YAML(typ='safe')
    for study in os.listdir(args.studies):
        try:
            with open(f'{args.studies}/{study}/{study}.yaml', 'r') as yaml_file:
                yaml = yaml_loader.load(yaml_file)
        except FileNotFoundError:
            print(f'{Fore.RED}{study}{Style.RESET_ALL}: No YAML file')
            continue
        except ruamel.yaml.parser.ParserError:
            print(f'{Fore.RED}{study}{Style.RESET_ALL}: Incorrectly formatted YAML file')

        # Special cases for flu
        if 'strain' in yaml and yaml['strain'] == 'Human adapted strain A/Aichi/2/1968, H3N2':
            species = 'H3N2_A_Aichi_2_1968'
        elif 'strain' in yaml and yaml['strain'] == 'Human adapted A/Perth/16/2009, H3N2':
            species = 'H3N2_A_Perth_16_2009'
        else:
            species = yaml['species'].lower().replace('. ', '_')

        name = f"{species}_{yaml['gene'].lower()}"

        if not name in seqs:
            print(f'{Fore.RED}{study}{Style.RESET_ALL}: No fasta sequence available')
        elif yaml['seq'] == seqs[name].seq:
            print(f'{Fore.GREEN}{study}{Style.RESET_ALL}: Sequence matches')
        else:
            print(f'{Fore.RED}{study}{Style.RESET_ALL}: Sequence does not match')


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--studies', '-s', help="Directory of study folders",
                        default=DEFAULT_STUDIES)

    parser.add_argument('--fasta', '-f', help="Directory of fasta files",
                        default=DEFAULT_FASTA)

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
