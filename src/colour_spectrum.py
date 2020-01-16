"""
Module containing functions to support creation of simple, continuous colour spectra
"""
import numpy as np
from  matplotlib import pyplot as plt
from matplotlib import cm

# TODO classes for types of colours?
# TODO way to determine output type to use
# TODO allow fewer values than the number of defined spectrum colours to be passes?

NA_COLOUR_DEFAULT = '0xC0C0C0'

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
    """
    def __init__(self, maximum, minimum, midpoint=None, colourmap=None, name='',
                 na_colour='0xC0C0C0'):
        self.maximum = maximum
        self.minimum = minimum

        self.midpoint = midpoint
        if self.midpoint is not None:
            diff = max(abs(self.maximum - self.midpoint), abs(self.minimum - self.midpoint))
            self.maximum = self.midpoint + diff
            self.minimum = self.midpoint - diff

        if colourmap is None:
            if midpoint is None:
                self.colourmap = cm.get_cmap('plasma')
            else:
                self.colourmap = cm.get_cmap('bwr')
        else:
            self.colourmap = cm.get_cmap(colourmap)

        self.na_colour = na_colour
        if not isinstance(na_colour, str):
            self.na_colour = rgb_to_hex(self.na_colour)

        self.name = name

    def __call__(self, val):
        rgb = self.colourmap(val - self.minimum) / (self.maximum - self.minimum)
        rgb = [rgb_clamp(x) for x in rgb[:3]]
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
