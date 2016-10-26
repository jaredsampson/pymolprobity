''' Utility functions for PyMOLProbity plugin. '''
import re


def quote_str(var, quote="'"):
    '''Enclose a str variable in quotes for printing.'''
    if type(var) is str:
        return '%s%s%s' % (quote, var, quote)
    else:
        return var


def slugify(s, sep='_'):
    '''Eliminate troublesome characters in a string.'''
    try:
        # Replace anything not alphanumeric with a separator.
        slug = re.sub(r'[^A-Za-z0-9]+', sep, s).strip(sep)

        # Condense multiple adjacent separators into one.
        mult_sep = r'[' + re.escape(sep) + ']+'
        slug = re.sub(mult_sep, sep, slug)

        return slug

    except TypeError:
        # Don't slugify non-string input.
        return s


def to_number(var):
    '''Convert a variable to a number if possible, and return it.

    Convert the passed variable to a number if possible and return the
    number.  Otherwise, return the original variable.
    '''
    try:
        # int (e.g. '1')
        return int(str(var))
    except:
        # float (e.g. '1.23')
        try:
            return float(str(var))
        except:
            # not a number
            return var

