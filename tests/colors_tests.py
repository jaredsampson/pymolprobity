import unittest

import mock

from .context import pymolprobity
import pymolprobity.colors as c


class GetPymolColorTests(unittest.TestCase):
    def test_with_color_list_color(self):
        inp = 'sea'
        ref = 'teal'
        res = c.get_pymol_color(inp)
        self.assertEqual(res, ref)

    def test_with_non_color_list_color(self):
        inp = 'somecolor'
        ref = 'somecolor'
        res = c.get_pymol_color(inp)
        self.assertEqual(res, ref)

    def test_with_None(self):
        inp = None
        ref = None
        res = c.get_pymol_color(inp)
        self.assertEqual(res, ref)


class GetColorRGBTests(unittest.TestCase):
    @mock.patch('pymolprobity.colors.cmd')
    def test_workflow(self, mock_cmd):
        mock_cmd.get_color_index.return_value = 0
        mock_cmd.get_color_tuple.return_value = (1.0, 1.0, 1.0)
        inp = 'white'
        ref = (1.0, 1.0, 1.0)
        res = c.get_color_rgb(inp)
        self.assertEqual(res, ref)
        mock_cmd.get_color_index.assert_called_with(inp)
        mock_cmd.get_color_tuple.assert_called_with(0)





