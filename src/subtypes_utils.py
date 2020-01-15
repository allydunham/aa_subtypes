#!/usr/bin/env python3
"""
Utility functions for the subtypes project
"""
from Bio.PDB.PDBIO import Select
from ruamel.yaml import YAML

def gene_to_filename(gene):
    """
    Convert a formal gene name to be suitable for file names
    """
    return gene.lower().replace(' ', '')

class SectionSelecter(Select):
    """
    Bio.PDBIO Selecter based on section input in the sytle of meta/structures.yaml
    """
    def __init__(self, sections, drop_hetero=False):
        self.sections = sections
        self.drop_hetero = drop_hetero

        for section in self.sections:
            if not 'region' in section:
                section['region'] = [float('-inf'), float('inf')]

    def accept_residue(self, residue):
        if self.drop_hetero and not residue.id[0] == ' ':
            return 0

        if self.sections is None:
            return 1

        chain = residue.full_id[2]
        position = residue.id[1]
        for section in self.sections:
            if (chain == section['chain'] and
                    section['region'][0] <= position <= section['region'][1]):
                return 1
        return 0

def import_sections(path_or_str, gene):
    """
    Import sections in the format of meta/structures.yaml. If path_or_string is a file it is
    assumed to be split into genes, each of which has a 'sections' field, whereas raw strings
    are assumed to be the sections of interest themselves, in YAML format.
    """
    yaml = YAML(typ='safe')

    if not path_or_str:
        return None

    try:
        with open(path_or_str, 'r') as yaml_file:
            sections = yaml.load(yaml_file)[gene]['sections']
    except FileNotFoundError:
        sections = yaml.load(path_or_str)

    return sections