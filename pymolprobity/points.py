'''Points (dots, vectors) for PyMOLProbity plugin.'''

import copy
import logging
import re

from chempy import cpv
#from pymol import cmd
from pymol import cgo

import colors
# from settings import mpgetq
import utils


logger = logging.getLogger(__name__)


###############################################################################
#
#  CGO UTILS
#
###############################################################################

def _cgo_color(color):
    '''Return a CGO list specifying a color.'''
    r, g, b = colors.get_color_rgb(color)
    return [cgo.COLOR, r, g, b]


def _cgo_sphere(pos, radius):
    '''Return a CGO list specifying a sphere.'''
    x, y, z = pos
    return [cgo.SPHERE, x, y, z, radius]


def _perp_vec(vec):
    '''Return a vector orthogonal to `vec`.'''
    if abs(vec[0]) > 1e-6:
        return [-(vec[1] + vec[2]) / vec[0], 1., 1.]
    if abs(vec[1]) > 1e-6:
        return [1., -(vec[0] + vec[2]) / vec[1], 1.]
    return [1., 1., -(vec[0] + vec[1]) / vec[2]]


def _cgo_quad(pos, normal, radius):
    '''Return a CGO list specifying a quad.'''
    v1 = cpv.normalize(_perp_vec(normal))
    v2 = cpv.cross_product(normal, v1)
    v1 = cpv.scale(v1, radius)
    v2 = cpv.scale(v2, radius)
    obj = [ cgo.BEGIN,
            cgo.TRIANGLE_STRIP,
            cgo.NORMAL]
    obj.extend(normal)
    obj.append(cgo.VERTEX)
    obj.extend(cpv.add(pos, v1))
    obj.append(cgo.VERTEX)
    obj.extend(cpv.add(pos, v2))
    obj.append(cgo.VERTEX)
    obj.extend(cpv.sub(pos, v2))
    obj.append(cgo.VERTEX)
    obj.extend(cpv.sub(pos, v1))
    obj.append(cgo.END)
    return obj


def _cgo_cylinder(pos0, pos1, radius, rgb0, rgb1):
    '''Return a CGO list specifying a cylinder.'''
    cgolist = [cgo.CYLINDER]
    cgolist.extend(pos0)    # start
    cgolist.extend(pos1)    # end
    cgolist.append(radius)
    cgolist.extend(rgb0)    # start color
    cgolist.extend(rgb1)    # end color
    return cgolist



##############################################################################
#
#  ATOM INFO STRINGS
#
###############################################################################

# KEY: N=name, A=alt, R=resn, I=resi, X=ins code, C=chain, O=occ, b=Bfac

# Simple dotlist atom: 15-char total, with 2-char chain ID
# e.g. "NNNNARRRIIIIXCC"
BASIC_ATOM_RE = re.compile(
        r"^([\w\? ]{4})"    # group 1: atom name (e.g. " CA ")
        r"([\w ])"          # group 2: alternate conformation ID (e.g. "A")
        r"([\w]{3})"        # group 3: residue name (e.g. "ARG")
        r"([\d ]{4}[\w ])"  # group 4: residue number + insertion code
        r"([\w ]{1,4})"     # group 5: chain ID (1-, 2-, or 4-char)
        )

def process_basic_atom_string(atom_str):
    m = BASIC_ATOM_RE.match(atom_str)
    if m:
        name, alt, resn, resi, chain = [g.strip().upper() for g in m.groups()]
        # Assemble into dict
        atom_dict = {
                'chain': chain,
                'resn': resn,
                'resi': resi,
                'name': name,
                'alt': alt,
                }
        return atom_dict
    else:
        return None


# Bonds vectorlist with occupancy (optional), B-factor, and input file name
# e.g. "NNNNARRRCCIIIIX OOOOBbbbbbb filename"
BONDS_VECTORLIST_ATOM_RE = re.compile(
        r"^([\w\? ]{4})"    # group 1: atom name (e.g. " CA ")
        r"([\w ])"          # group 2: alternate conformation ID (e.g. "A")
        r"([\w]{3})"        # group 3: residue name (e.g. "ARG")
        r"([\w ]{2})"       # group 4: chain ID
        r"([\d ]{4}[\w ]) " # group 5: residue number + ins code
        r"([\d\.]{4}|) ?"   # group 6: 4-char occupancy (optional) + space?
        r"B([\d\.]{4,6}) "  # group 7: 4- to 6-char B-factor + space
        r".+$"              # input file name
        )

def process_bonds_vectorlist_atom(atom_str):
    m = BONDS_VECTORLIST_ATOM_RE.match(atom_str)
    if m:
        (name, alt, resn, chain,
                resi, occ, b) = [g.strip().upper() for g in m.groups()]
        # Convert occupancy and B-factor to floats
        if not occ:
            occ = 1.00
        else:
            occ = float(occ)
        if not b:
            b = 0.00
        else:
            b = float(b)
        # Assemble into dict
        atom_dict = {
                'chain': chain,
                'resn': resn,
                'resi': resi,
                'name': name,
                'alt': alt,
                'occ': occ,
                'b': b,
                }
        return atom_dict
    else:
        return None



###############################################################################
#
#  DOTS
#
###############################################################################

class Dot(object):
    """Python representation of a Molprobity-style kinemage dot."""

#     def set_draw(self):
#         self.draw = 1

#     def unset_draw(self):
#         self.draw = 0

#     def toggle_draw(self):
#         if self.draw:
#             self.draw = 0
#         else:
#             self.draw = 1

    def get_cgo(self, dot_mode=0, dot_radius=0.03):
        """Generate a CGO list for a dot."""
        cgolist = []

        # COLOR
        cgolist.extend(_cgo_color(self.color))

        if dot_mode == 0:  # spheres
            logger.debug("Adding dot to cgolist...")
            cgolist.extend(_cgo_sphere(self.coords, dot_radius))
            logger.debug("Finished adding dot to cgolist.")

        if dot_mode == 1:  # quads
            logger.debug("Adding quad to cgolist...")
            normal = cpv.normalize(cpv.sub(self.coords, self.atom['coords']))
            cgolist.extend(_cgo_quad(self.coords, normal, dot_radius * 1.5))
            logger.debug("Finished adding quad to cgolist.")

        return cgolist

    def __init__(self,
            atom=None, color=None, pointmaster=None, coords=None, draw=1):
        # Atom
        self.atom = atom

        # Dot info
        self.color = colors.get_pymol_color(color)
        self.pm = pointmaster
        self.coords = coords
        self.draw = draw

        # Dotlist info
        self.dotlist_name = None
        self.dotlist_color = None
        self.master = None

#         # Bind to PyMOL atom
#         #self.atom_selection = None
#         #self.atom_id = None
#         #if self.atom is not None:
#             #self.atom_selection = self._get_atom_selection()
#             #self.atom_id = cmd.id_atom(self.atom_selection)
#             #print self.atom_selection, self.atom_id


DOTLIST_HEADER_RE = re.compile(
        r"dotlist "         # dotlist keyword + space
        r"{([^}]*)} "       # atom info section + space
        r"color=([^\s]*) "  # color
        r"master={([^}]*)}" # master
        )

def _parse_dotlist_header(line):
    """Parse the header line of a kinemage `@dotlist` keyword.

    Header lines are in the following format:

        dotlist {x} color=white master={vdw contact}
        dotlist {x} color=sky master={vdw contact}
        dotlist {x} color=red master={H-bonds}

    Where "x" is an arbitrary name for the dotlist (currently hard-coded as "x"
    in the Probe source), the color is the default color of the dots, and
    master is the type of interaction depicted by the dots.

    """
    m = DOTLIST_HEADER_RE.match(line)

    # name, color, master
    return  m.group(1), m.group(2), utils.slugify(m.group(3))


DOTLIST_BODY_RE = re.compile(
            r"{([^}]*)}"        # atom info string
            r"(\w*)\s*"         # color + optional space(s)
            r"'(\w)' "          # pointmaster
            r"([0-9.,\-]*)"     # coordinates
            )

def _parse_dotlist_body(lines):
    """Parse the non-header lines of a kinemage `@dotlist` keyword.

    Body lines are in the following format[*]:

        { CA  SER  26  A}blue  'O' 61.716,59.833,8.961
        { OE1 GLU  31  A}blue  'S' 61.865,58.936,17.234
        { OE2 GLU  31  A}blue  'S' 61.108,60.399,15.044
        { H?  HOH 293  A}greentint  'O' 57.884,59.181,7.525 [**]
        { O   HOH 435  A}greentint  'O' 56.838,61.938,21.538
        { O   HOH 450  A}greentint  'O' 55.912,56.611,17.956

    Or generally:

        {AAAABCCCDDDDEFF}colorname  'G' X,Y,Z

    where:

        AAAA = atom name
        B    = alt conf
        CCC  = residue name (3-letter)
        DDDD = residue number (typically 3 digits)
        E    = insertion code
        FF   = chain (typically a single letter)
        G    = pointmaster code (e.g. ScSc, McSc, McMc, Hets)

    * Note that the text within the braces is not space-delimited, but is
    arranged in fixed-width columns.

    ** Also, 'H?' for water hydrogens is problematic and simply stripped away.

    """
    active_atom = None
    dots = []

    # Parse lines to generate Dots
    for i, l in enumerate(lines):
        m = DOTLIST_BODY_RE.match(l)

        # Atom selection in kinemages is only written explicitly the first
        # time for a given set of dots. Afterward, it inherits from the
        # previous line via a single double-quote character (").
        if m.group(1) == '"':
            atom = active_atom
        else:
            # e.g.: " O   HOH 450  A"
            logger.debug('m.group(1) match is: %s' % m.group(1))
            atom_sel = m.group(1)
            # TODO: Check this formatting in probe documentation
            name = atom_sel[0:4].strip().replace('?','')
            alt = atom_sel[4:5].strip().upper()
            resn = atom_sel[5:8].strip().upper()
            resi = atom_sel[8:13].strip().upper()
            chain = atom_sel[13:15].strip().upper()

            logger.debug('Generating atom for dot %i...' % i)
            # TODO don't create duplicate atoms (track in MPObject)
            atom = {'name': name,
                    'alt': alt,
                    'resn': resn,
                    'resi': resi,
                    'chain': chain}
            logger.debug('Finished generating atom for dot %i.' % i)

            active_atom = atom

        color = m.group(2)
        pointmaster = m.group(3)
        coords = [float(c) for c in m.group(4).split(',')]

        # Create the Dot
        dot = Dot(atom, color, pointmaster, coords)
        dots.append(dot)

    return dots


# #def _get_dot_atom_selection(dot):
#     #assert type(dot) is Dot
#     #obj = "%s" % dot.dotlist.result.obj
#     #if dot.atom['chain']:
#         #sele = "%s and chain %s" % (sele, dot.chain)
#     #if dot.atom['resn']:
#         #sele = "%s and resn %s" % (sele, dot.resn)
#     #if dot.atom['resi']:
#         #sele = "%s and resi %s" % (sele, dot.resi)

#     ## Hack: Probe gives HOH hydrogens atom names of 'H?', which, even when the
#     ## '?' is stripped, doesn't work with PyMOL, which numbers them 'H1' and
#     ## 'H2'.
#     #if dot.atom['name'] and not dot.atom['resn'] == 'HOH':  # hack
#         #sele = "%s and name %s" % (sele, dot.name)
#     #return sele


def process_dotlist(lines, context):
    '''Process a list of dotlist lines and return a list of Dots.

    Given a list of lines from a Kinemage file comprising a dotlist, parse the
    first line as the header, and the remaining lines as the body.  Create a
    Dot() instance for each line and return a list of Dots.

    '''
    logger.debug("Parsing dotlist header...")
    name, color, master = _parse_dotlist_header(lines[0])
    logger.debug("Parsing dotlist body...")
    dots = _parse_dotlist_body(lines[1:])

    logger.debug("Adding Dotlist info...")
    # Add Dotlist info to each dot.
    for d in dots:
        # From dotlist header
        d.dotlist_name = name
        d.dotlist_color = color
        d.master = master

        # From context
        d.kinemage = context['kinemage']
        d.group = context['group']
        d.subgroup = context['subgroup']
        d.animate = context['animate']

    logger.debug("Finished adding Dotlist info.")

    return dots




# ###############################################################################
# #
# #  VECTORS
# #
# ###############################################################################

class Vector(object):
    """Python representation of a Molprobity-style kinemage vector."""

#     def set_draw(self):
#         self.draw = 1

#     def unset_draw(self):
#         self.draw = 0

#     def toggle_draw(self):
#         if self.draw:
#             self.draw = 0
#         else:
#             self.draw = 1

    def get_cgo(self, radius=0.03):
        """Generate a CGO list for a vector."""
        cgolist = []

        # Set colors
        rgb0 = colors.get_color_rgb(self.color[0])
        rgb1 = colors.get_color_rgb(self.color[1])

        if True:  # cylinders
            logger.debug("Adding vector to cgolist...")
            # Cylinder
            cgolist.extend(_cgo_cylinder(self.coords[0], self.coords[1],
                radius, rgb0, rgb1))
            # Caps
            cgolist.extend(_cgo_color(self.color[0]))
            cgolist.extend(_cgo_sphere(self.coords[0], radius))
            cgolist.extend(_cgo_color(self.color[1]))
            cgolist.extend(_cgo_sphere(self.coords[1], radius))

            logger.debug("Finished adding vector to cgolist.")

        return cgolist

    def macro(self, i):
        a = copy.copy(self.atom[i])
        if a is None:
            return None
        if a['alt']:
            a['alt'] = '`{}'.format(a['alt'])
        else:
            a['alt'] = ''
        return '{chain}/{resn}`{resi}/{name}{alt}'.format(**a)

    def sel(self, i):
        a = copy.copy(self.atom[i])
        if a is None:
            return None
        c = 'chain {chain} and '.format(**a) if a['chain'] else ''
        i = 'resi {resi} and '.format(**a) if a['resi'] else ''
        n = 'name {name} and '.format(**a) if a['name'] else ''
        alt = 'alt {alt} and '.format(**a) if a['alt'] else ''

        sel = '{}{}{}{}'.format(c, i, n, alt)
        return sel[:-5]  # strip last " and "


    def __init__(self,
            atom0=None, color0=None, pointmaster0=None, coords0=None,
            atom1=None, color1=None, pointmaster1=None, coords1=None, draw=1):
        # Atom
        self.atom = [atom0, atom1]

        # Vector info
        c0 = colors.get_pymol_color(color0)
        c1 = colors.get_pymol_color(color1)
        self.color = [c0, c1]
        self.pm = [pointmaster0, pointmaster1]
        self.coords = [coords0, coords1]
        self.draw = draw

        # Vectorlist info
        self.vectorlist_name = None
        self.vectorlist_color = None
        self.master = None

#         # Bind to PyMOL atom
#         #self.atom_selection = None
#         #self.atom_id = None
#         #if self.atom is not None:
#             #self.atom_selection = self._get_atom_selection()
#             #self.atom_id = cmd.id_atom(self.atom_selection)
#             #print self.atom_selection, self.atom_id

    def __str__(self):
        vectorlist_info = '[{},{}]'.format(self.vectorlist_name, self.master)
        return '{}: {}--{}'.format(vectorlist_info, self.macro(0), self.macro(1))


VECTORLIST_HEADER_RE = re.compile(
        r"vectorlist "
        r"{([^}]*)} "           # group 1: name
        r"color= *([^\s]*) *"   # group 2: color
        r"(\w+ )*"              # group 3: other text (e.g. nobutton)
        r"master= *{([^}]*)}"   # group 4: master
        )
def _parse_vectorlist_header(line):
    """Parse the header line of a kinemage `@vectorlist` keyword.

    Header lines are in the following format:
        vectorlist {x} color=white master={small overlap}

    Where "x" is an arbitrary name for the list (currently hard-coded as "x"
    in the Probe source), the color is the default color of the vectors, and
    master is the type of interaction depicted by the vectors.

    """
    logger.debug('parsing line: "{}"'.format(line))

    NAME = 1
    COLOR = 2
    OTHER = 3
    MASTER = 4
    m = VECTORLIST_HEADER_RE.match(line)

    logger.debug('vectorlist header: {}'.format(m.groups()))

    # name, color, master
    return m.group(NAME), m.group(COLOR), utils.slugify(m.group(MASTER))


# Probe v2.16 (20-May-13) format
VECTORLIST_CLASH_RE = re.compile(
        r"{([^}]*)}"        # group 1: atom description
        r"(\w*) "           # group 2: color
        r"([A-Z] )*"        # group 3: optional L or P character followed by space
        r" *"               # allow extra space before pointmaster
        r"'(\w)' "          # group 4: pointmaster
        r"([0-9.,\-]*)"     # group 5: coordinates as a single string
        )
def _parse_clash_vectorlist_body(lines):
    """Parse the non-header lines of a vectorlist containing clash spikes.

    Body lines are in the following format:
        { CB  SER  26  A}yellowtint P  'O' 57.581,59.168,8.642 {"}yellowtint   'O' 57.589,59.163,8.646

    Or generally:

        {AAAA BBB CCCD E}colorname  'F' X,Y,Z       (x 2)

    where:

        AAAA = atom name
        BBB  = residue name (3-letter)
        CCC  = residue number
        D    = insertion code  ??? TODO: check this
        E    = chain
        F    = pointmaster code (e.g. ScSc, McSc, McMc, Hets)

    * Note that the text within the braces is not space-delimited, but is
    arranged in fixed-width columns.

    Also, 'H?' for water hydrogens is problematic.  # TODO

    """

    ATOM = 1
    COLOR = 2
    LP = 3
    PM = 4
    COORDS = 5

    active_atom = None
    vectors = []

    # Parse lines to generate Vectors
    for i, l in enumerate(lines):
        matches = VECTORLIST_CLASH_RE.finditer(l)

        v = []
        for m in matches:
            logger.debug('match: {}'.format(m.group(0)))
            logger.debug('clash vectorlist body line: {}'.format(m.groups()))
            logger.debug('beginning match...')
            # Atom selection in kinemages is only written explicitly the first
            # time for a given list of points. Afterward, it inherits from the
            # previous point via a single double-quote character (").
            atom_sel = m.group(ATOM)
            if atom_sel == '"':
                #logger.debug('using active atom...')
                atom = active_atom
            else:
                # TODO don't create duplicate atoms (track in MPObject)
                logger.debug('Generating atom for vector %i point...' % i)
                atom = process_basic_atom_string(atom_sel)
                logger.debug('Finished generating atom for vector %i point.' % i)

                active_atom = atom

            color = m.group(COLOR)
            pointmaster = m.group(PM)
            coords = [float(c) for c in m.group(COORDS).split(',')]

            v.append({
                    'atom': atom,
                    'color': color,
                    'pm': pointmaster,
                    'coords': coords})

        # Create the Vector
        vector = Vector(v[0]['atom'], v[0]['color'], v[0]['pm'], v[0]['coords'],
                        v[1]['atom'], v[1]['color'], v[1]['pm'], v[1]['coords'])
        vectors.append(vector)

    return vectors


VECTORLIST_BONDS_RE = re.compile(
        r"{([^}]*)}"        # group 1: atom description
        r" ?"               # optional space
        r"([LP]) "          # group 2: single L or P character (TODO: what is this?)
        r"(?:'(\w)' )*"     # group 3: pointmaster (optional)
        r"(?:(\w+) )*"      # group 4: other non-quoted word, e.g. 'ghost' (optional)
        r"([\d,\.\- ]+)"    # group 5: coordinates as a single string
        )

def _parse_bonds_vectorlist_body(lines):
    '''Parse vectorlist lines that describe bonds.'''
    # TODO: merge this with _parse_clashes_vectorlist_body()

    ATOM = 1
    LP = 2      # not used
    PM = 3
    OTHER = 4   # not used
    COORDS = 5

    # active_atom = None
    vectors = []

    for i, l in enumerate(lines):
        matches = tuple(VECTORLIST_BONDS_RE.finditer(l))

        # Set up new vector points list
        v = []

        # If a continuation point, reuse the second point from the previous vector
        if len(matches) == 1:
            try:
                prev_v = vectors[-1]
                v.append({
                        'atom': prev_v.atom[1],
                        'color': prev_v.color[1],
                        'pm': prev_v.pm[1],
                        'coords': prev_v.coords[1]})
            except IndexError:
                logger.error('only 1 point given, but no previous points to use')
                raise
        elif len(matches) < 1 or len(matches) > 2:
            # <1 or >2 shouldn't happen
            raise ValueError('malformed line: {}'.format(l))

        for j, m in enumerate(matches):
            logger.debug('match {} of {}: {}'.format(j+1, len(matches), m.group(0)))
            logger.debug('bond vectorlist body line: {}'.format(m.groups()))
            logger.debug('beginning match...')
            # # Atom selection in kinemages is only written explicitly the first
            # # time for a given set of dots. Afterward, it inherits from the
            # # previous line via a single double-quote character (").
            atom_sel = m.group(ATOM)
            if atom_sel == '"':
                logger.debug('using active atom...')
                atom = active_atom
            else:
                try:
                    # TODO don't create duplicate atoms (track in MPObject)
                    msg = 'Generating atom for vector {} point {}...'.format(i, j)
                    logger.debug(msg)

                    atom = process_bonds_vectorlist_atom(atom_sel)
                except:
                    msg = 'Atom info string `{}` could not be parsed. Skipping.'
                    logger.error(msg.format(atom_sel))

            logger.debug('Finished generating atom for vector %i point %i.' % (i, j))

            # active_atom = atom

            color = None
            pointmaster = m.group(PM)
            coords = [float(c) for c in m.group(COORDS).split(',')]

            v.append({
                    'atom': atom,
                    'color': color,
                    'pm': pointmaster,
                    'coords': coords})

        # logger.debug('coords0: {}'.format(v[0]['coords']))
        # logger.debug('coords1: {}'.format(v[1]['coords']))

        # Create the Vector
        vector = Vector(v[0]['atom'], v[0]['color'], v[0]['pm'], v[0]['coords'],
                        v[1]['atom'], v[1]['color'], v[1]['pm'], v[1]['coords'])

        vectors.append(vector)

    return vectors




# #def _get_dot_atom_selection(dot):
#     #assert type(dot) is Dot
#     #obj = "%s" % dot.dotlist.result.obj
#     #if dot.atom['chain']:
#         #sele = "%s and chain %s" % (sele, dot.chain)
#     #if dot.atom['resn']:
#         #sele = "%s and resn %s" % (sele, dot.resn)
#     #if dot.atom['resi']:
#         #sele = "%s and resi %s" % (sele, dot.resi)

#     ## Hack: Probe gives HOH hydrogens atom names of 'H?', which, even when the
#     ## '?' is stripped, doesn't work with PyMOL, which numbers them 'H1' and
#     ## 'H2'.
#     #if dot.atom['name'] and not dot.atom['resn'] == 'HOH':  # hack
#         #sele = "%s and name %s" % (sele, dot.name)
#     #return sele


def process_vectorlist(lines, context):
    """Process a list of dotlist lines and return a list of Vectors.

    Given a list of lines from a Kinemage file comprising a dotlist, parse the
    first line as the header, and the remaining lines as the body.  Create a
    Dot() instance for each line and return a list of Dots.

    """
    logger.debug("Parsing vectorlist header...")
    name, color, master = _parse_vectorlist_header(lines[0])
    logger.debug("Parsing vectorlist body...")
    if name == 'x':
        vectors = _parse_clash_vectorlist_body(lines[1:])
    else:
        vectors = _parse_bonds_vectorlist_body(lines[1:])

    logger.debug("Adding vectorlist info...")
    # Add Dotlist info to each dot.
    for v in vectors:
        # From vectorlist header
        v.vectorlist_name = name
        v.vectorlist_color = color
        v.master = master

        # From context
        v.kinemage = context['kinemage']
        v.group = context['group']
        v.subgroup = context['subgroup']
        v.animate = context['animate']

    logger.debug("Finished adding Vectorlist info.")

    return vectors


