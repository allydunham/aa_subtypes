"""
Module containing functions to support creation of simple, continuous colour spectra
"""
from bisect import bisect

# TODO classes for types of colours?
# TODO way to determine output type to use
# TODO allow fewer values than the number of defined spectrum colours to be passes?

NA_COLOUR_DEFAULT = '0xC0C0C0'
SPECTRA = {
    'RdBu': [(178, 24, 43), (214, 96, 77), (244, 165, 130), (253, 219, 199),
             (247, 247, 247), (209, 229, 240), (146, 197, 222), (67, 147, 195),
             (33, 102, 172)]
}

class ColourSpectrum:
    """
    Define a colour spectrum, which returns the appropriate Hex codes when called with
    a numeric value. The base constructor provides a versatible but manual definition system,
    but most uses are better suited to one of the alternate constructors (e.g. cls.linear()).

    values:    Iterable of numeric values giving the position of each colour
    colours:   String spectrum name (as given in colour_spectrum.SPECTRA) or an iterable
               of RGB tuples giving the colour at each value.
    na_colour: RGB tuple or Hex code of colour to return for out of range values
    """
    def __init__(self, values, colours='RdBu', na_colour='0xC0C0C0'):
        try:
            # If colours doesn't index SPECTRA assume its an iterable of colours
            spectrum = colours
            colours = SPECTRA[colours]
        except TypeError:
            spectrum = None

        self.spectrum = spectrum

        self.na_colour = na_colour
        if not isinstance(na_colour, str):
            self.na_colour = rgb_to_hex(self.na_colour)

        self.max_value = max(values)
        self.min_value = min(values)

        if len(values) != len(colours):
            raise ValueError(("values and colours must be the same length."
                              "This means the appropriate number of values must be provided "
                              f"when using a built-in scale. In this case {len(colours)} colours "
                              "were given"))

        # Sort values and colours so interpolation makes sense
        values = sorted(list(zip(values, colours)), key=lambda x: x[0])
        self.values, self.colours = zip(*values)

    @classmethod
    def linear(cls, minimum, maximum, midpoint=None, colours='RdBu', na_colour=NA_COLOUR_DEFAULT):
        """
        Create a linear spectrum interpolating between bounds, with optional symmetry about
        a midpoint.

        maximum: Upper limit of the scale
        minimum: Lower limit of the scale
        midpoint: Midpoint to make the scale symmetric about, causing maximum or minimum to
                  be adjusted such that they are equidistant from midpoint
        colour: String name of a spectrum (as given in colour_spectrum.SPECTRA) or an iterable of
                RGB tuples to interpolate between
        na_colour: RGB tuple or Hex code of colour to return for out of range values
        """
        if not maximum > minimum:
            raise ValueError('maximum must be greater than minimum')

        if colours in SPECTRA.keys():
            colours = SPECTRA[colours]

        n_colours = len(colours)

        # Use midpoint to automatically adjust max and min
        if midpoint is not None:
            if n_colours % 2 == 0:
                mid_colour = rgb_interpolate(colours[n_colours // 2 - 1],
                                             colours[n_colours // 2],
                                             0.5)
                colours.insert(n_colours // 2, mid_colour)

            diff = max(abs(maximum - midpoint), abs(minimum - midpoint))
            maximum = midpoint + diff
            minimum = midpoint - diff

        colour_div = (maximum - minimum)/(n_colours - 1)
        values = [minimum + colour_div * x for x in range(n_colours)]
        return cls(values, colours, na_colour)

    def __call__(self, val):
        if val in self.values:
            return rgb_to_hex(self.colours[self.values.index(val)])

        ind = bisect(self.values, val)

        if ind == 0 or ind >= len(self.values):
            return self.na_colour

        prop = (val - self.values[ind - 1]) / (self.values[ind] - self.values[ind - 1])
        res = rgb_interpolate(self.colours[ind - 1], self.colours[ind], prop)

        return rgb_to_hex(res)

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
