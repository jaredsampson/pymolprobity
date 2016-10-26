'''Tests for PyMOLProbity flips module.'''

import mock
import unittest

from .context import pymolprobity
import pymolprobity.flips as flips



class TestFlip(unittest.TestCase):
    def test_base_init(self):
        res = flips.Flip()
        self.assertEqual(res.flip_class, None)
        self.assertEqual(res.set, None)
        self.assertEqual(res.set_index, None)
        self.assertEqual(res.chain, None)
        self.assertEqual(res.resi, None)
        self.assertEqual(res.resn, None)
        self.assertEqual(res.name, None)
        self.assertEqual(res.alt, None)
        self.assertEqual(res.flip_type, None)
        self.assertEqual(res.descr, None)
        self.assertEqual(res.best_score, None)
        self.assertEqual(res.best_has_bad_bump, None)
        self.assertEqual(res.init_score, None)
        self.assertEqual(res.init_has_bad_bump, None)
        self.assertEqual(res.code, None)
        self.assertEqual(res.orig_score, None)
        self.assertEqual(res.orig_has_bad_bump, None)
        self.assertEqual(res.flip_score, None)
        self.assertEqual(res.flip_has_bad_bump, None)
        self.assertEqual(res.reduce_flipped, None)
        self.assertEqual(res.user_flipped, None)

    def test_macro(self):
        f = flips.Flip()
        f.chain = 'A'
        f.resi = '100'
        f.resn = 'THR'
        f.name = 'CA'
        res = f.macro()
        ref = 'A/THR`100/CA'
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.flips.Flip.macro')
    def test_str(self, mock_macro):
        mock_macro.return_value = 'macro'
        f = flips.Flip()
        res = f.__str__()
        ref = 'Flip: macro'
        self.assertEqual(res, ref)


@mock.patch('pymolprobity.flips.parse_other_score')
@mock.patch('pymolprobity.flips.parse_canflip_score')
@mock.patch('pymolprobity.flips.parse_methyl_score')
@mock.patch('pymolprobity.flips.parse_func_group')
class TestParseFlips(unittest.TestCase):
    def setUp(self):
        self.set = ['Set10.1: B 358 ASN     :      amide:sc=   -2.71! C(o=-1.3!,f=-3.4!)']
        self.single = ['Single : A  41 ASN     :FLIP  amide:sc= -0.0583  X(o=-0.064,f=-0.058)']

    def test_user_mod_with_too_few_sections(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        user_mod = ['Set:2']
        with self.assertRaises(ValueError):
            flips.parse_flips(user_mod)

    def test_user_mod_with_too_many_sections(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        user_mod = ['Set:2:3:4:5']
        with self.assertRaises(ValueError):
            flips.parse_flips(user_mod)

    def test_with_set(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_canflip_score.return_value = [1, 2, 3, 4, 5, 6, 7]
        user_mod = self.set
        res = flips.parse_flips(user_mod)
        assert type(res) is list
        assert len(res) == 1
        f = res[0]
        assert f.flip_class == 'set'
        assert f.set == '10'
        assert f.set_index == '1'

    def test_with_single(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_canflip_score.return_value = [1, 2, 3, 4, 5, 6, 7]
        user_mod = self.single
        res = flips.parse_flips(user_mod)
        assert type(res) is list
        assert len(res) == 1
        f = res[0]
        assert f.flip_class == 'single'
        assert f.set is None
        assert f.set_index is None

    def test_with_non_single_non_set(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        '''Parsing USER MOD with no flip lines should return an empty list.'''
        user_mod = ['blah']  # not Set or Single
        res = flips.parse_flips(user_mod)
        mock_pfunc.assert_not_called()
        mock_parse_methyl_score.assert_not_called()
        mock_parse_canflip_score.assert_not_called()
        mock_parse_other_score.assert_not_called()
        self.assertEqual(res, [])

    def test_with_set_not_matching_set_number_regex(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        user_mod = ['Set blah:2:3:4']
        with self.assertRaises(ValueError):
            flips.parse_flips(user_mod)


    def test_calls_parse_func_group(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_canflip_score.return_value = [1, 2, 3, 4, 5, 6, 7]
        user_mod = self.single
        res = flips.parse_flips(user_mod)
        mock_pfunc.assert_called_once_with(' A  41 ASN     ')

    def test_sets_flip_func_group_info(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = ['A', '100B', 'THR', 'OG1', '']
        mock_parse_canflip_score.return_value = [1, 2, 3, 4, 5, 6, 7]
        user_mod = self.single
        res = flips.parse_flips(user_mod)
        f = res[0]
        assert f.chain == 'A'
        assert f.resi == '100B'
        assert f.resn == 'THR'
        assert f.name == 'OG1'
        assert f.alt == ''

    def test_single_methyl_calls_parse_methyl_score(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_methyl_score.return_value = [1, 2, 3, 4]
        user_mod = ['Single : A  57 LYS NZ  :NH3+    180:sc=       0   (180deg=0)']
        flips.parse_flips(user_mod)
        assert mock_parse_methyl_score.called

    def test_single_canflip_calls_parse_canflip_score(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_canflip_score.return_value = [1, 2, 3, 4, 5, 6, 7]
        user_mod = ['Single : A  41 ASN     :FLIP  amide:sc= -0.0583  X(o=-0.064,f=-0.058)']
        flips.parse_flips(user_mod)
        assert mock_parse_canflip_score.called

    def test_single_other_calls_parse_other_score(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        mock_parse_other_score.return_value = [1, 2]
        user_mod = ['Single : A  89 SER OG  :   rot  180:sc=       0']
        flips.parse_flips(user_mod)
        assert mock_parse_other_score.called

    def test_unparseable_score(self, mock_pfunc,
            mock_parse_methyl_score, mock_parse_canflip_score,
            mock_parse_other_score):
        mock_pfunc.return_value = [1, 2, 3, 4, 5]
        user_mod = ['Single : A  89 SER OG  :   rot  180:bad score section']
        with self.assertRaises(ValueError):
            flips.parse_flips(user_mod)


class TestParseFuncGroup(unittest.TestCase):
    def test_with_1_char_chain_ids(self):
        fg = 'B 510 HIS     '  # len 14
        ref = ['B','510','HIS','','']
        res = flips.parse_func_group(fg)
        self.assertEqual(res, ref)

    def test_with_2_char_chain_ids(self):
        fg = ' B 510 HIS     '  # len 15
        ref = ['B','510','HIS','','']
        res = flips.parse_func_group(fg)
        self.assertEqual(res, ref)

    def test_with_4_char_chain_ids(self):
        fg = '   B 510 HIS     '  # len 17
        ref = ['B','510','HIS','','']
        res = flips.parse_func_group(fg)
        self.assertEqual(res, ref)

    def test_with_res_and_name(self):
        fg = ' B 528 THR OG1 '
        ref = ['B', '528', 'THR', 'OG1', '']
        res = flips.parse_func_group(fg)
        self.assertEqual(res, ref)

    def test_with_altconf(self):
        fg = ' B 528 THR OG1A'
        ref = ['B', '528', 'THR', 'OG1', 'A']
        res = flips.parse_func_group(fg)
        self.assertEqual(res, ref)

    def test_with_bad_input_length(self):
        with self.assertRaises(ValueError):
            flips.parse_func_group('blah')


class ParseMethylTests(unittest.TestCase):
    def test_with_basic_input(self):
        groups = ('       0', ' ', '0', '')
        ref = (0.0, False, 0.0, False)
        res = flips.parse_methyl_score(groups)
        self.assertEqual(res, ref)

    def test_with_actual_numbers_and_1_clash(self):
        groups = ('   0.728', ' ', '-1.25', '!')
        ref = (0.728, False, -1.25, True)
        res = flips.parse_methyl_score(groups)
        self.assertEqual(res, ref)


class ParseCanflipTests(unittest.TestCase):
    def test_with_basic_input(self):
        groups = (' -0.0583', ' ', 'X', '-0.064', '', '-0.058', '')
        ref = (-0.0583, False, 'X', -0.064, False, -0.058, False)
        res = flips.parse_canflip_score(groups)
        self.assertEqual(res, ref)

    def test_with_some_clashes(self):
        groups = ('   -2.71', '!', 'C', '-1.3', '!', '-3.4', '!')
        ref = (-2.71, True, 'C', -1.3, True, -3.4, True)
        res = flips.parse_canflip_score(groups)
        self.assertEqual(res, ref)


class ParseOtherTests(unittest.TestCase):
    def test_with_basic_input(self):
        groups = ('       0', '')
        ref = (0, False)
        res = flips.parse_other_score(groups)
        self.assertEqual(res, ref)

    def test_with_number_and_clash(self):
        groups = ('    -1.2', '!')
        ref = (-1.2, True)
        res = flips.parse_other_score(groups)
        self.assertEqual(res, ref)




if __name__ == '__main__':
    unittest.main()

