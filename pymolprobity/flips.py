'''Flips for PyMOLProbity plugin.'''

from __future__ import print_function

import logging
import re

# from pymol import cmd

logger = logging.getLogger(__name__)


class Flip(object):
    """Store information related to a bond rotation optimized by Reduce."""
    def macro(self):
        '''Return an atom selection macro for the Flip.'''
        c = self.chain or ''
        i = self.resi or ''
        n = self.name or ''
        r = self.resn or ''
        sep = '`' if (i and r) else ''
        macro = '{}/{}{}{}/{}'.format(c, r, sep, i, n)
        return macro

    def __init__(self):
        self.flip_class = None
        self.set = None
        self.set_index = None

        self.chain = None
        self.resi = None
        self.resn = None
        self.name = None
        self.alt = None

        self.flip_type = None
        self.descr = None
        self.best_score = None
        self.best_has_bad_bump = None

        # for "methyl" type (H added to freely rotatable bond)
        self.init_score = None
        self.init_has_bad_bump = None

        # for "canflip" type (e.g. NQH flips)
        self.code = None
        self.orig_score = None
        self.orig_has_bad_bump = None
        self.flip_score = None
        self.flip_has_bad_bump = None

        # Track user choices for making/reverting flips
        self.reduce_flipped = None
        self.user_flipped = None

    def __str__(self):
        return 'Flip: {}'.format(self.macro())


# Based on: http://stackoverflow.com/a/4914089/501277
def slices(input_str, col_lengths):
    '''Split a string into columns of the lengths in the col_lengths list.'''
    position = 0
    for length in col_lengths:
        yield input_str[position:position + length]
        position += length


def parse_func_group(func_group):
    '''Parse functional group info from second section of a flip.

    Returns 5 strings: chain, resi, resn, name, alt

    Functional group strings are written by Reduce via the source code:
    ```
    ::sprintf(descrbuf, "%-2.2s%4s%c%-3.3s%-4.4s%c",
        hr.chain(), hr.Hy36resno(), hr.insCode(),
        hr.resname(), c1.atomname(), hr.alt());
    std::string descr = descrbuf;
    ```

    and generally look something like, e.g. " A 254 THR OG1 ".

    '''
    # Set columns based on above source code excerpt, with residue number and
    # insertion code merged into the 2nd column
    col_lengths = [2, 5, 3, 4, 1]

    # Default is 15-character string.
    if len(func_group) != 15:
        # Allow for different chain ID lengths (which can apparently happen,
        # according to the `flipkin` script).
        if len(func_group) == 14:
            # 1-character chain ID field
            col_lengths[0] = 1
        elif len(func_group) == 17:
            # 4-character chain ID field
            col_lengths[0] = 4
        else:
            # Non-standard functional group string format
            msg = 'Non-standard length ({}) of functional group string: `{}`'
            raise ValueError(msg.format(len(func_group), func_group))

    # Slice it up & strip whitespace
    return [col.strip() for col in list(slices(func_group, col_lengths))]


def parse_methyl_score(groups):
    '''Parse the score section of a methyl flip match.'''
    best_score = float(groups[0].strip())
    best_has_bad_bump = (groups[1] == '!')
    init_score = float(groups[2])
    init_has_bad_bump = (groups[3] == '!')
    return best_score, best_has_bad_bump, init_score, init_has_bad_bump


def parse_canflip_score(groups):
    '''Parse the score section of a canflip flip match.'''
    best_score = float(groups[0])
    best_has_bad_bump = (groups[1] == '!')
    code = groups[2]
    orig_score = float(groups[3])
    orig_has_bad_bump = (groups[4] == '!')
    flip_score = float(groups[5])
    flip_has_bad_bump = (groups[6] == '!')
    return (best_score, best_has_bad_bump, code, orig_score,
            orig_has_bad_bump, flip_score, flip_has_bad_bump)


def parse_other_score(groups):
    '''Parse the score section of an other flip match.'''
    best_score = float(groups[0])
    best_has_bad_bump = (groups[1] == '!')
    return best_score, best_has_bad_bump


def parse_flips(user_mod):
    '''Parse an array of USER MOD lines and return a list of Flips.'''
    ret = []
    set_re = re.compile(r"Set\s?([0-9]*)\.([0-9])")
    for line in user_mod:
        if line.startswith('Single') or line.startswith('Set'):
            try:
                flip_class, func_group, descr, score = line.split(':')
            except ValueError:
                logger.error("Malformed USER MOD line: `{}`".format(line))
                raise
        else:  # pragma: nocover
            continue

        # new flip
        f = Flip()

        # Determine flip class
        if flip_class.startswith('Single'):
            f.flip_class = 'single'
        else:
            set_match = set_re.match(flip_class)
            if set_match:
                f.flip_class = 'set'
                g = set_match.groups()
                f.set = g[0]
                f.set_index = g[1]
            else:
                msg = 'Malformed Set flip: {}'.format(flip_class)
                raise ValueError(msg)

        # Parse functional group (i.e. residue/atom info)
        f.chain, f.resi, f.resn, f.name, f.alt = parse_func_group(func_group)

        # Store orientation description as-is
        f.descr = descr

        # Determine flip type (canflip, methyl, other) from score section
        methyl_score = re.compile(
                r"sc=(.{8})"              # best score
                r"(!| )  "  # 2 spaces    # bad clash indicator for best score
                r"\(180deg=([\-\.0-9]*)"  # initial score
                r"(!?)\)"                 # bad clash indicator (may be empty)
                )
        canflip_score = re.compile(
                r"sc=(.{8})"              # best score
                r"(!| ) "   # 1 space     # bad clash indicator for best score
                r"(.)"                    # class code (Reduce's recommendation)
                r"\(o=([\-\.0-9]*)"       # original max score
                r"(!?),"                  # bad clash indicator for original
                r"f=([\-\.0-9]*)"         # flipped max score
                r"(!?)\)"                 # bad clash indicator for flipped
                )
        other_score = re.compile(
                # Note: This matches the preceding two types of flip records,
                # and therefore must be checked last.
                r"sc=(.{8})"              # best score
                r"(!?)"                   # optional bad clash indicator
                )

        # Match to regex, then process the score section.
        # Methyl
        m = methyl_score.match(score)
        if m:
            f.flip_type = 'methyl'
            (f.best_score, f.best_score_has_bad_bump, f.init_score,
                    f.init_has_bad_bump) = parse_methyl_score(m.groups())
        else:
            # Canflip
            m = canflip_score.match(score)
            if m:
                f.flip_type = 'canflip'
                (f.best_score, f.best_score_has_bad_bump,
                    f.code, f.orig_score, f.orig_has_bad_bump, f.flip_score,
                    f.flip_has_bad_bump) = parse_canflip_score(m.groups())
                f.reduce_flipped = 1 if (f.code == 'F') else 0
            else:
                # Other
                m = other_score.match(score)
                if m:
                    f.flip_type = 'other'
                    (f.best_score,
                        f.best_score_has_bad_bump) = parse_other_score(m.groups())
                else:
                    raise ValueError("UNABLE TO PARSE FLIP!\n"
                                 "raw flip record:\n"
                                 ">>{}<<".format(user_mod))
        ret.append(f)

    return ret
