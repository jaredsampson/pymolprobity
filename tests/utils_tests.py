'''utils tests for pymolprobity plugin'''

import mock
import unittest

from .context import pymolprobity
import pymolprobity.utils as ut


class QuoteStrTests(unittest.TestCase):
    def test_with_string_input(self):
        res = ut.quote_str("inp")
        ref = "'inp'"
        self.assertEqual(res, ref)

    def test_with_nonstring_input(self):
        res = ut.quote_str(1)
        self.assertEqual(res, 1)

    def test_with_custom_quote(self):
        res = ut.quote_str("inp", '?')
        ref = "?inp?"
        self.assertEqual(res, ref)


class SlugifyTests(unittest.TestCase):
    def test_basic(self):
        inp = 'Some string 1'
        res = ut.slugify(inp)
        ref = 'Some_string_1'
        self.assertEqual(res, ref)

    def test_leading_trailing_junk(self):
        inp = ' this string had leading and trailing space '
        res = ut.slugify(inp)
        ref = 'this_string_had_leading_and_trailing_space'
        self.assertEqual(res, ref)

    def test_multiple_adjacent_separators(self):
        inp = 'this   string   had   extra    spaces'
        res = ut.slugify(inp)
        ref = 'this_string_had_extra_spaces'
        self.assertEqual(res, ref)

    def test_custom_separator(self):
        inp = 'some string'
        res = ut.slugify(inp, sep='!')
        ref = 'some!string'
        self.assertEqual(res, ref)

    def test_non_string_input(self):
        inp = 1
        res = ut.slugify(inp)
        ref = 1
        self.assertEqual(res, ref)


class ToNumberTests(unittest.TestCase):
    def test_with_int_str(self):
        inp = '1'
        res = ut.to_number(inp)
        ref = 1
        self.assertEqual(res, ref)

    def test_with_float_str(self):
        inp = '1.23'
        res = ut.to_number(inp)
        ref = 1.23
        self.assertEqual(res, ref)

    def test_with_scientific_notation_str(self):
        inp = '1e-2'
        res = ut.to_number(inp)
        ref = 0.01
        self.assertEqual(res, ref)

    def test_with_non_numeric_str(self):
        inp = 'foo'
        res = ut.to_number(inp)
        ref = inp
        self.assertEqual(res, ref)



if __name__ == '__main__':
    unittest.main()

