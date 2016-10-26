'''Color functions for PyMOLProbity plugin.'''

from pymol import cmd


def get_pymol_color(color):
    """Return the PyMOL color corresponding to a Kinemage color name."""

    color_list = {
        # Color names defined in the KiNG source that aren't included in a
        # standard PyMOL installation are listed here with the names of
        # (approximately) equivalent PyMOL named colors.

        #"king_color": "pymol_color",
        "sea":         "teal",
        "sky":         "skyblue",
        "peach":       "yelloworange",
        "lilac":       "arsenic",
        "pinktint":    "lightpink",
        "peachtint":   "lightorange",
        "yellowtint":  "paleyellow",
        "greentint":   "palegreen",
        "bluetint":    "lightblue",
        "lilactint":   "lithium",
        "deadwhite":   "white",
        "deadblack":   "black",
        }

    try:
        return color_list[color]
    except:
        return color


def get_color_rgb(color):
    """Given a color name, returns the RGB values."""
    index = cmd.get_color_index(color)
    rgb = cmd.get_color_tuple(index)
    return rgb  # e.g. (1.0, 1.0, 1.0)

