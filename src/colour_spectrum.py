"""
Module containing functions to support creation of simple, colour spectra with continuous
and discrete options. Currently setup to produce PyMol style hexcodes (0xRRGGBB)
"""
import re
import numpy as np
from  matplotlib import pyplot as plt
from matplotlib import cm

NA_COLOUR_DEFAULT = '0xDCDCDC'

RE_HEXCODE = re.compile('[A-Za-z0-9]{6}$')

class ColourPalette:
    """
    Simple map between categories and colours, with the same interface as ColourSpectrum.
    Input colours can

    colourmap: Dictionary of value: colour pairs
    name:      Pallete name
    na_colour: Colour to return when not in colourmap
    """
    def __init__(self, colourmap, name='', na_colour=NA_COLOUR_DEFAULT, order=None):
        self.name = name
        self.na_colour = make_hex(na_colour)

        self.colourmap = {k: make_hex(v) for k, v in colourmap.items()}

        self.order = order or sorted(self.colourmap.keys())

    def __call__(self, value):
        return self.colourmap.get(value, self.na_colour)

    def plot(self, horizontal=False):
        """
        Plot the palette as a free legend
        """
        colours = [hex_to_rgb(self.colourmap[i]) for i in self.order]
        height = [1 for _ in self.order]
        size = (2, 0.5) if horizontal else (0.5, 2)

        fig, axis = plt.subplots(figsize=size)
        axis.set_frame_on(False)
        axis.get_yaxis().set_visible(False)
        axis.get_xaxis().set_visible(False)
        axis.tick_params(axis='both', which='both', bottom=False,
                         top=False, left=False, right=False)

        if horizontal:
            axis.set_xlabel(self.name)
            axis.scatter(x=self.order, y=height, c=colours)
        else:
            axis.set_ylabel(self.name)
            axis.get_yaxis().set_label_position("right")
            axis.scatter(x=height, y=self.order, c=colours)

        return fig, axis

class ColourSpectrum:
    """
    Wrap a matplotlib colourmap as a callable returning Hexcodes for given values, for
    use in manual plotting (for instance colouring PyMOL structures)

    maximum:   Upper limit of the scale
    minimum:   Lower limit of the scale
    midpoint:  Midpoint to make the scale symmetric about, causing maximum or minimum to
               be adjusted such that they are equidistant from midpoint
    colourmap: Matplotlib colourmap defining the main spectrum
    na_colour: RGB tuple or Hex code of colour to return for out of range values
    na_outside_range: Use na_colour for values not between minimum or maximum, otherwise
                      defer this to the colourmap.
    """
    def __init__(self, minimum, maximum, midpoint=None, colourmap=None, name='',
                 na_colour=NA_COLOUR_DEFAULT, na_outside_range=True):
        self.maximum = maximum
        self.minimum = minimum

        self.midpoint = midpoint
        if self.midpoint is not None:
            diff = max(abs(self.maximum - self.midpoint), abs(self.midpoint - self.minimum))
            self.maximum = self.midpoint + diff
            self.minimum = self.midpoint - diff

        if colourmap is None:
            if midpoint is None:
                self.colourmap = cm.get_cmap('plasma')
            else:
                self.colourmap = cm.get_cmap('bwr')
        else:
            self.colourmap = cm.get_cmap(colourmap)

        self.na_colour = make_hex(na_colour)
        self.na_outside_range = na_outside_range

        self.name = name

    def __call__(self, val):
        if np.isnan(val):
            return self.na_colour

        if self.na_outside_range and not self.minimum <= val <= self.maximum:
            return self.na_colour

        rgb = self.colourmap((val  - self.minimum) / (self.maximum - self.minimum))
        rgb = [rgb_clamp(x * 255) for x in rgb[:3]]
        return rgb_to_hex(rgb)

    def plot(self, horizontal=False):
        """
        Plot a colourbar of the spectrum as a free plot
        """
        div = (self.maximum - self.minimum) / 100
        image = np.arange(self.minimum, self.maximum + div, div)
        if horizontal:
            image = np.vstack([image, image])
            size = (2, 0.5)
            extent = (self.minimum, self.maximum, 0, 1)
        else:
            image = np.flip(image)
            image = np.vstack([image, image]).T
            size = (0.5, 2)
            extent = (0, 1, self.minimum, self.maximum)

        fig, axis = plt.subplots(figsize=size)
        axis.imshow(image, aspect='auto', cmap=self.colourmap, extent=extent)
        axis.set_frame_on(False)

        if horizontal:
            axis.get_yaxis().set_visible(False)
            axis.set_xlabel(self.name)

        else:
            axis.get_xaxis().set_visible(False)
            axis.get_yaxis().tick_right()
            axis.get_yaxis().set_label_position("right")
            axis.set_ylabel(self.name)

        return fig, axis

def rgb_interpolate(low, high, prop):
    """
    Interpolate between two RGB tuples
    """
    return [rgb_clamp(y + (x - y) * prop) for x, y in zip(low, high)]

def rgb_clamp(colour_value):
    """
    Clamp a value to integers on the RGB 0-255 range
    """
    return int(min(255, max(0, colour_value)))

def rgb_to_hex(rgb):
    """
    Covert RGB tuple to Hex code
    """
    return f'0x{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}'

def hex_to_rgb(hexcode):
    """
    Convert Hex code to RGB tuple
    """
    return (int(hexcode[-6:-4], 16), int(hexcode[-4:-2], 16), int(hexcode[-2:], 16))

def make_hex(value):
    """
    Convert to PyMOl style hexcodes
    """
    if isinstance(value, str):
        if RE_HEXCODE.search(value):
            value = f'0x{value.upper()[-6:]}'
        else:
            raise ValueError(f'Unrecognised string number: {value}. Enter a hexcode or rgb triplet')

    elif isinstance(value, (tuple, list)) and len(value) == 3:
        if all(0 < x < 1 for x in value):
            value = [rgb_clamp(x * 255) for x in value]
        value = rgb_to_hex(value)

    else:
        raise ValueError(f'Unrecognised number format: {value}. Enter a hexcode or rgb triplet')

    return value
