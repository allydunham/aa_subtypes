#!/usr/bin/env python3
"""
Utility functions for the subtypes project
"""

from io import StringIO

import pandas as pd
import numpy as np
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB import PDBParser
from Bio.Alphabet.IUPAC import IUPACProtein
from ruamel.yaml import YAML

from colour_spectrum import ColourSpectrum, ColourPalette

AMINO_ACID_COLOURS = {
    'A': '0xFF0000', 'I': '0xFA8072', 'L': '0xB22222', 'M': '0xFF8C00', 'V': '0xFF6347',
    'F':'0xFFD700', 'W': '0xyellow3', 'Y': '0xDAA520', 'N': '0xB8860B', 'C': '0x6495ED',
    'Q': '0x00FFFF', 'S': '0x0000FF', 'T': '0x483D8B', 'R': '0x00FF00', 'H': '0x228B22',
    'K': '0x98FB98', 'D': '0xBA55D3', 'E': '0xFFC0CB', 'G': '0xFAEBD7',
    'P': '0x000000', 'X': '0x808080'
}

class SubtypesColourMap:
    """
    Colour mappings for deep mutational landscape propertues, looking up colours
    and adjusting scales specifically for this projects factors
    """
    # Discrete palettes. (name, colourmap, na colour)
    palettes = {
        'wt': ('WT AA', AMINO_ACID_COLOURS, '0xDCDCDC', None),
        'mut': ('Mutant AA', AMINO_ACID_COLOURS, '0xDCDCDC', None),
        'cluster_num': ('Subtype', {
            0: '0xDCDCDC', 1: '0xE41A1C', 2: '0x377EB8', 3: '0x4DAF4A',
            4: '0x984EA3', 5: '0xFF7F00', 6: '0xFFFF33', 7: '0xA65628',
            8: '0xF781BF'
        }, '0x808080', None)
    }

    # (name, scale, midpoint, na_colour) for properties
    # default is ($property, 'bwr', 0, '0xDCDCDC')
    spectra = {
        'PC1': ('PC1', 'coolwarm', 0, '0xDCDCDC'),
        'PC2': ('PC2', 'BrBG', 0, '0xDCDCDC'),
        'PC3': ('PC3', 'PuOr', 0, '0xDCDCDC'),
        'PC4': ('PC4', 'RdYlGn', 0, '0xDCDCDC'),
        'tSNE1': ('tSNE1', 'viridis', 0, '0x29AF7F'),
        'tSNE2': ('tSNE2', 'plasma', 0, '0xDE7065'),
        'umap1': ('UMAP1', 'viridis', 0, '0x29AF7F'),
        'umap2': ('UMAP2', 'plasma', 0, '0xDE7065'),
        'total_energy': ('Mean FoldX Energy', 'PiYG', 0, '0xDCDCDC'),
        'mean_sift': ('Mean log10 SIFT', 'PuRd_r', None, '0xDCDCDC'),
        'mean_score': (f'Mean Norm. ER', 'RdBu', 0, '0xDCDCDC')
    }
    spectra = dict(**spectra,
                   **{aa: (f'Norm. ER ({aa})', 'RdBu', 0, '0xDCDCDC') for
                      aa in IUPACProtein.letters})

    # Exclude 1% outliers from scale (give same as max colour)
    outlier_properties = ['PC1'] + list(IUPACProtein.letters)

    # Explicitally clamp some properties with known interpretations (e.g. energy from FoldX)
    clamped_properties = {'total_energy': (-7.5, 7.5)}

    @classmethod
    def lookup_map(cls, landscape_property, values=None):
        """
        Fetch the appropriate map for a given property
        """
        if landscape_property in cls.palettes.keys():
            return SubtypesColourPalette(landscape_property)

        else:
            return SubtypesColourSpectrum(landscape_property, values)

class SubtypesColourPalette(ColourPalette, SubtypesColourMap):
    """
    Special colour palettes for deep mutational landscape properties
    """
    def __init__(self, landscape_property):
        name, colourmap, na_colour, order = self.palettes[landscape_property]
        super().__init__(colourmap=colourmap, name=name, na_colour=na_colour, order=order)

class SubtypesColourSpectrum(ColourSpectrum, SubtypesColourMap):
    """
    Special colour spectra for deep mutational landscape properties, looking up spectra
    and adjusting scales specifically for this projects factors
    """
    def __init__(self, landscape_property, values):
        na_outside = True
        values = np.array(values)
        values = values[~np.isnan(values)]
        minimum = min(values)
        maximum = max(values)

        spec = self.spectra.get(landscape_property, (landscape_property, 'bwr', 0, '0xDCDCDC'))
        name, spectrum, midpoint, na_colour = spec

        if landscape_property in self.outlier_properties:
            minimum = np.quantile(values, 0.01)
            maximum = np.quantile(values, 0.99)
            na_outside = False

        elif landscape_property in self.clamped_properties:
            clamp_min, clamp_max = self.clamped_properties[landscape_property]
            minimum = max(minimum, clamp_min)
            maximum = min(maximum, clamp_max)
            na_outside = False

        super().__init__(minimum, maximum, midpoint=midpoint, colourmap=spectrum, name=name,
                         na_colour=na_colour, na_outside_range=na_outside)

def gene_to_filename(gene):
    """
    Convert a formal gene name to be suitable for file names
    """
    return gene.lower().replace(' ', '')

def import_sections(path_or_str, gene):
    """
    Import sections in the format of meta/structures.yaml. If path_or_string is a file it is
    assumed to be split into genes, each of which has a 'sections' field, whereas raw strings
    are assumed to be the sections of interest themselves, in YAML format.
    """
    yaml = YAML(typ='safe')

    try:
        with open(path_or_str, 'r') as yaml_file:
            sections = yaml.load(yaml_file)[gene]['sections']
    except FileNotFoundError:
        sections = yaml.load(path_or_str)

    return sections

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

def get_pdb_string(pdb_path, select=None):
    """
    Get string representation of a PDB file, filtered using select as in Bio.PDB.PDBIO
    """
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('_', pdb_path)

    pdbio = PDBIO()
    pdbio.set_structure(structure)
    with StringIO('PDB') as virtual_pdb_file:
        pdbio.save(virtual_pdb_file, select=select)
        return virtual_pdb_file.getvalue()

def offset_uniprot_position(position, sections):
    """
    Convert Uniprot positions to PDB positions given a list of sections determining the
    possible offsets. Each section dict must have a chain, integer offset (PDB -> uniprot) and
    two item region item.
    """
    for sec in sections:
        reg = sec['region'] if 'region' in sec else [float('-inf'), float('inf')]
        if reg[0] <= position - sec['offset'] <= reg[1]:
            return sec['chain'], position - sec['offset']

    return np.nan, np.nan

def dms_pdb_positions(dms, sections):
    """
    Get PDB positions for all positions in a Deep Scanning Dataframe, with columns for gene,
    and uniprot position, based on a dict of sections for each gene (can also be in the
    structures.yaml format)
    """
    try:
        # Try to convert structures.yaml format to section dict
        sections = {k: v['sections'] for k, v in sections.items()}
    except TypeError:
        # Otherwise assume it is already a dict of section lists
        pass

    chains = []
    pdb_positions = []

    for gene, position in zip(dms.gene, dms.position):
        cha, pos = offset_uniprot_position(position, sections[gene_to_filename(gene)])
        chains.append(cha)
        pdb_positions.append(pos)

    return chains, pdb_positions

def quick_load(cluster=None):
    """
    Helper function to load sections, dms, etc. for interactive session
    """
    yaml = YAML(typ='safe')
    with open('meta/structures.yaml', 'r') as yaml_file:
        sections = {k: v['sections'] for k, v in yaml.load(yaml_file).items()}

    dms = pd.read_csv('data/combined_mutational_scans.tsv', sep='\t')
    dms['pdb_chain'], dms['pdb_position'] = dms_pdb_positions(dms, sections)

    if cluster is not None:
        cluster = pd.read_csv(cluster, sep='\t')
        dms = pd.merge(cluster, dms, how='left', on=['study', 'gene', 'position', 'wt'])
        dms['cluster_num'] = dms['cluster'].slice(1).astype(int)

    return sections, dms
