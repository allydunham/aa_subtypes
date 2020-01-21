"""
PyMOL commands to aid producing figures of the mutational landscape projected
onto proteins.
"""
from pymol import cmd
from colour_spectrum import ColourSpectrum

def colour_spectrum(chains, positions, values, colourer=None, na_colour='0xC0C0C0'):
    """Colour specific residues according to a colourmap. colourer must return a Hexcode
    when called with a value as well as have an 'na_colour' attribute if no na_colour
    is specifically supplied.
    """
    if colourer is None:
        colourer = ColourSpectrum(min(values), max(values), colourmap='viridis')

    cmd.color(na_colour or colourer.na_colour, 'prot')

    for chn, pos, val in zip(chains, positions, values):
        cmd.color(colourer(val), f'prot and chain {chn} and resi {int(pos)}')
