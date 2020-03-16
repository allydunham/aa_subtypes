"""
Library for visualising data in PyMOL, paticularly projecting arbitary values onto proteins.
"""
from itertools import cycle

from colour_spectrum import ColourSpectrum
import subtypes_utils as su

def quick_highlight(cmd, dms, gene, factor='cluster_num', res=None, root='data/pdb'):
    """
    Load a new protein and highlight a given factor.
    """
    cmd.delete('all')
    cmd.load(f'{root}/{gene}.pdb')

    if res:
        pdb_dms = dms[(dms.gene == gene) & (dms.wt == res)]
    else:
        pdb_dms = dms[dms.gene == gene]
    pdb_dms = pdb_dms.dropna(subset=['pdb_position'])

    colourer = su.SubtypesColourMap.lookup_map(factor, dms[factor])
    project_landscape(cmd, pdb_dms.pdb_chain, pdb_dms.pdb_position, pdb_dms.cluster_num, colourer)
    return pdb_dms

def project_landscape(cmd, chain, position, value, colourer=None, na_colour=None):
    """
    Colour specific residues according to a colourmap. colourer must return a Hexcode
    when called with a value as well as have an 'na_colour' attribute if no na_colour
    is specifically supplied. Chain can either be a single identifier (str) or an
    iterable of identifiers
    """
    if colourer is None:
        colourer = ColourSpectrum(min(value), max(value), colourmap='viridis')

    if isinstance(chain, str):
        chain = cycle([chain])

    colour_residues(cmd, *zip(chain, position, [colourer(val) for val in value]),
                    base_colour=na_colour or colourer.na_colour)

def colour_residues(cmd, *args, base_colour=None):
    """
    Colour multiple residues programatically. Each argument should be a
    (chain, position index, hex code) tuple
    """
    if base_colour is not None:
        cmd.color(base_colour, 'polymer')

    for chn, pos, col in args:
        pos = int(pos)
        pos = f'\{pos}' if pos < 0 else pos # Negative indices must be escaped
        cmd.color(col, f'polymer and chain {chn} and resi {pos}')
