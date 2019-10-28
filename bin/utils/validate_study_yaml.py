#!/usr/bin/env python3
"""
Check study config YAML files (assumed to be in {study}/{study}.yaml format) for expected format,
consistent sequences (Fasta files are assumed to be named in the format
{species e.g. h_sapiens}_{gene}.fa/fasta except for the special case of influenza) etc.
"""
# TODO Rework messages to be much more standardised, maybe a class based system?
# Look for existing implementations

import argparse
import os
from collections import defaultdict

import ruamel.yaml
from Bio import SeqIO
from colorama import Fore
from colorama import Style

import subtypes_utils as sutil

PROJECT_ROOT = '/Users/ally/phd/subtypes'
DEFAULT_STUDIES = f'{PROJECT_ROOT}/data/studies'
DEFAULT_FASTA = f'{PROJECT_ROOT}/meta/fasta'

EXPECTED_FIELDS = ('study', 'gene', 'uniprot_id', 'gene_type', 'species', 'seq',
                   'experiment', 'transform', 'authour', 'year', 'title', 'doi',
                   'pmid', 'url')

HEADERS = ['Study', 'YAML Exists?', 'ID Correct?', 'Seq Correct?', 'Missing Fields?']

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

def check_yaml_seq(yaml, seqs):
    """
    Check a config YAML seq against a list of reference seqs
    Returns: Match, No fasta or Mismatch
    """
    # Get species
    if yaml['species'] == 'Influenza':
        if yaml['strain'] == 'Human adapted strain A/Aichi/2/1968, H3N2':
            species = 'H3N2_A_Aichi_2_1968'
        elif yaml['strain'] == 'Human adapted A/Perth/16/2009, H3N2':
            species = 'H3N2_A_Perth_16_2009'
    else:
        species = yaml['species'].lower().replace('. ', '_')

    name = f"{species}_{sutil.gene_to_filename(yaml['gene'])}"

    if not name in seqs:
        return 'No fasta'

    elif yaml['seq'] == seqs[name].seq:
        return 'Match'

    else:
        return 'Mismatch'

def check_study(study, root_dir='.', seqs=None):
    """
    Validate a study. Returns a dictionary containing the validation data:
    yaml: Yes or an error message
    missing_fields: a list of missing fields
    id: Yes or an error message
    seq: Match, No fasta or Mismatch
    """
    yaml_loader = ruamel.yaml.YAML(typ='safe')
    output = defaultdict(str)
    # Check YAML
    try:
        with open(f'{root_dir}/{study}/{study}.yaml', 'r') as yaml_file:
            yaml = yaml_loader.load(yaml_file)
        output['yaml'] = 'Yes'

    except FileNotFoundError:
        output['yaml'] = 'No YAML'
        return output

    except ruamel.yaml.parser.ParserError:
        output['yaml'] = 'Bad YAML'
        return output

    # Check ID
    output['id'] = 'Yes' if yaml['study'] == study else 'Mismatch'

    # Check Seq
    if seqs is not None:
        output['seq'] = check_yaml_seq(yaml, seqs)

    # Check fields
    output['missing_fields'] = [x for x in EXPECTED_FIELDS if not x in yaml.keys()]

    return output

def print_study_validation_table(studies):
    """
    Print a formated validation table from a dictionary of check_study output
    """
    try:
        print(Style.BRIGHT, Fore.BLUE, sep='', end='')
        print("{: <25} {: <15} {: <15} {: <15} {: <20}".format(*HEADERS))
        for study, checks in studies.items():
            colours = defaultdict(lambda: Fore.GREEN)

            if not checks['yaml'] == 'Yes':
                colours['yaml'] = Fore.RED
                colours['study'] = Fore.RED

            if not checks['id'] == 'Yes':
                colours['id'] = Fore.RED
                colours['study'] = Fore.RED

            if not checks['seq'] == 'Match':
                colours['seq'] = Fore.RED
                colours['study'] = Fore.RED

            if checks['missing_fields'] == []:
                missing_fields = 'None'
            else:
                missing_fields = ', '.join(checks['missing_fields'])
                colours['missing_fields'] = Fore.RED
                colours['study'] = Fore.RED

            output = [f"{colours['study']}{study}",
                      f"{colours['yaml']}{checks['yaml']}",
                      f"{colours['id']}{checks['id']}",
                      f"{colours['seq']}{checks['seq']}",
                      f"{colours['missing_fields']}{missing_fields}"]
            print("{: <30} {: <20} {: <20} {: <20} {: <20}".format(*output))
        print(Style.RESET_ALL)
    except:
        print(Style.RESET_ALL)
        raise

def main(args):
    """Main script"""
    seqs = import_single_seq_fasta_dir(args.fasta)
    studies = {s: check_study(s, args.studies, seqs) for s in os.listdir(args.studies)}
    print_study_validation_table(studies)

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
