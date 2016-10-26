import copy
import unittest

import mock

from .context import pymolprobity
import pymolprobity.kinemage as kin



class KinemageTests(unittest.TestCase):
    def setUp(self):
        self.kw_types = ['viewid', 'group', 'subgroup', 'master',
                'pointmaster', 'dotlist', 'vectorlist']
        self.kin = kin.Kinemage()
        for i, kw in enumerate(self.kw_types):
            self.kin.keywords[i] = {'keyword': kw,
                                    'data': '{} data'.format(kw)}
        # duplicate each item to test get_unique
        self.kin2 = copy.deepcopy(self.kin)
        for i in range(0, 7):
            self.kin2.keywords[i-7] = self.kin2.keywords[i]

    def tearDown(self):
        self.kin = None

    def test_get_all_keywords_of_type(self):
        for kw in self.kw_types:
            res = self.kin.get_all_keywords_of_type(kw)
            self.assertEqual(len(res), 1)
            self.assertEqual(res[0], '{} data'.format(kw))

    def test_get_unique_keywords_of_type(self):
        for kw in self.kw_types:
            res = self.kin2.get_unique_keywords_of_type(kw)
            self.assertEqual(len(res), 1)
            self.assertEqual(res[0], '{} data'.format(kw))

    @mock.patch('pymolprobity.kinemage.Kinemage.get_all_keywords_of_type')
    def test_viewids(self, mock_get_all):
        res = self.kin.viewids()
        ref = mock_get_all.return_value
        mock_get_all.assert_called_once_with('viewid')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_all_keywords_of_type')
    def test_groups(self, mock_get_all):
        res = self.kin.groups()
        ref = mock_get_all.return_value
        mock_get_all.assert_called_once_with('group')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_unique_keywords_of_type')
    def test_subgroups(self, mock_get_unique):
        res = self.kin.subgroups()
        ref = mock_get_unique.return_value
        mock_get_unique.assert_called_once_with('subgroup')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_unique_keywords_of_type')
    def test_masters(self, mock_get_unique):
        res = self.kin.masters()
        ref = mock_get_unique.return_value
        mock_get_unique.assert_called_once_with('master')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_unique_keywords_of_type')
    def test_pointmasters(self, mock_get_unique):
        res = self.kin.pointmasters()
        ref = mock_get_unique.return_value
        mock_get_unique.assert_called_once_with('pointmaster')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_all_keywords_of_type')
    def test_dotlists(self, mock_get_all):
        res = self.kin.dotlists()
        ref = mock_get_all.return_value
        mock_get_all.assert_called_once_with('dotlist')
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.kinemage.Kinemage.get_all_keywords_of_type')
    def test_vectorlists(self, mock_get_all):
        res = self.kin.vectorlists()
        ref = mock_get_all.return_value
        mock_get_all.assert_called_once_with('vectorlist')
        self.assertEqual(res, ref)


class KinemageDrawMethodTests(unittest.TestCase):
    # TODO
    pass


class ProcessKinemageTests(unittest.TestCase):
    def setUp(self):
        self.context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    @mock.patch('pymolprobity.kinemage.points.process_dotlist')
    def test_calls_process_dotlist_with_dotlists(self,
            mock_proc_dotlist):
        inp = '@dotlist blah'
        mock_proc_dotlist.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_dotlist.assert_called_once_with(
                ['dotlist blah'], self.context)
        ref = {'keyword': 'dotlist', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.points.process_vectorlist')
    def test_calls_process_vectorlist_with_vectorlists(self,
            mock_proc_vectorlist):
        inp = '@vectorlist blah'
        mock_proc_vectorlist.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_vectorlist.assert_called_once_with(
                ['vectorlist blah'], self.context)
        ref = {'keyword': 'vectorlist', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_viewid')
    def test_calls_process_viewid_with_viewids(self,
            mock_proc_viewid):
        inp = '@viewid blah'
        mock_proc_viewid.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_viewid.assert_called_once_with(
                ['viewid blah'], self.context)
        ref = {'keyword': 'viewid', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_master')
    def test_calls_process_master_with_master(self,
            mock_proc_master):
        inp = '@master blah'
        mock_proc_master.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(
                ['master blah'], self.context)
        ref = {'keyword': 'master', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_pointmaster')
    def test_calls_process_pointmaster_with_pointmaster(self,
            mock_proc_pm):
        inp = '@pointmaster blah'
        mock_proc_pm.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_pm.assert_called_once_with(
                ['pointmaster blah'], self.context)
        ref = {'keyword': 'pointmaster', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_kinemage_keyword')
    def test_calls_process_kinemage_keyword_with_kinemage(self,
            mock_proc_kin):
        inp = '@kinemage blah'
        mock_proc_kin.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc_kin.assert_called_once_with(
                ['kinemage blah'], self.context)
        ref = {'keyword': 'kinemage', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_group')
    def test_calls_process_group_with_group(self,
            mock_proc):
        inp = '@group blah'
        mock_proc.return_value = 'val'
        k = kin.process_kinemage(inp)
        mock_proc.assert_called_once_with(
                ['group blah'], self.context)
        ref = {'keyword': 'group', 'data': 'val'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.process_subgroup')
    def test_calls_process_subgroup_with_subgroup(self,
            mock_proc):
        inp = '@subgroup blah'
        context = {}
        mock_proc.return_value = 'data'
        k = kin.process_kinemage(inp)
        mock_proc.assert_called_once_with(
                ['subgroup blah'], self.context)
        ref = {'keyword': 'subgroup', 'data': 'data'}
        self.assertEqual(k.keywords[0], ref)

    @mock.patch('pymolprobity.kinemage.points.process_vectorlist')
    @mock.patch('pymolprobity.kinemage.points.process_dotlist')
    def test_with_skipped_keyword(self, mock_proc_dotlist,
            mock_proc_vectorlist):
        inp = '@text something'
        kin.process_kinemage(inp)
        mock_proc_dotlist.assert_not_called()
        mock_proc_vectorlist.assert_not_called()
        # TODO: test prints debug message

    @mock.patch('pymolprobity.kinemage.logger')
    def test_with_unknown_keyword(self, mock_logger):
        inp = '@not_a_keyword blah'
        k = kin.process_kinemage(inp)
        mock_logger.warning.assert_called_with('Unknown keyword: not_a_keyword')

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_kinemage_keyword')
    def test_kinemage_keyword_updates_context(self, mock_proc_kin,
            mock_proc_master):
        inp = '@kinemage blah\n@master blah'
        context = self.context
        context['kinemage'] = mock_proc_kin.return_value
        kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(['master blah'], context)

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_group')
    def test_group_updates_context(self, mock_proc_group,
            mock_proc_master):
        inp = '@group blah\n@master blah'
        context = self.context
        context['group'] = mock_proc_group.return_value
        kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(['master blah'], context)

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_group')
    def test_none_group_updates_context(self, mock_proc_group,
            mock_proc_master):
        inp = '@group blah\n@group blah\n@master blah'
        mock_proc_group.side_effect = ( ['reduce', 'animate'], None )

        context1 = copy.deepcopy(self.context)
        context1['group'] = ['reduce', 'animate']
        context1['animate'] = 1

        context2 = copy.deepcopy(context1)
        context2['group'] = None  # from 2nd group
        context2['animate'] = 0

        kin.process_kinemage(inp)

        mock_proc_group.assert_has_calls(
                [mock.call(['group blah'], self.context),
                 mock.call(['group blah'], context1)])
        mock_proc_master.assert_called_once_with(['master blah'], context2)

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_group')
    def test_animate_group_updates_context(self, mock_proc_group,
            mock_proc_master):
        inp = '@group blah\n@master blah'
        mock_proc_group.return_value = ['reduce', 'animate']
        context = self.context
        context['group'] = mock_proc_group.return_value
        context['animate'] = 1
        kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(['master blah'], context)

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_group')
    def test_non_animate_group_updates_context(self, mock_proc_group,
            mock_proc_master):
        inp = '@group blah\n@group blah\n@master blah'
        # first call sets animate = 1, second should reset it to 0
        mock_proc_group.side_effect = [['animate'], ['blah']]
        context = self.context
        context['group'] = ['blah']
        context['animate'] = 0
        kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(['master blah'], context)

    @mock.patch('pymolprobity.kinemage.process_master')
    @mock.patch('pymolprobity.kinemage.process_subgroup')
    def test_subgroup_updates_context(self, mock_proc_subgroup,
            mock_proc_master):
        inp = '@subgroup blah\n@master blah'
        context = self.context
        context['subgroup'] = mock_proc_subgroup.return_value
        kin.process_kinemage(inp)
        mock_proc_master.assert_called_once_with(['master blah'], context)






@mock.patch('pymolprobity.kinemage.logger')
class SingleLineKeywordCheckTests(unittest.TestCase):
    def test_with_single_line(self, mock_logger):
        inp = ['line 1']
        kin.single_line_keyword_check(inp)
        self.assertFalse(mock_logger.warning.called)

    def test_with_multiple_lines(self, mock_logger):
        inp = ['line 1', 'line 2']
        kin.single_line_keyword_check(inp)
        self.assertTrue(mock_logger.warning.called)

    def test_with_non_list_input(self, mock_logger):
        inp = 42
        with self.assertRaises(ValueError):
            kin.single_line_keyword_check(inp)


class ProcessViewidTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    @mock.patch('pymolprobity.kinemage.single_line_keyword_check')
    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_viewid(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_first_viewid(self):
        inp = ['viewid { Q28   A}']
        res = kin.process_viewid(inp, self.base_context)
        ref = {
                'view_num': 1,
                'flipped': False,
                'resn': 'Q',
                'resi': '28',
                'alt': '',
                'chain': 'A',
                }
        self.assertEqual(res, ref)

    def test_with_second_or_later_viewid(self):
        inp = ['2viewid { Q32   A}']
        res = kin.process_viewid(inp, self.base_context)
        self.assertEqual(res['view_num'], 2)

    def test_with_3_digit_resnum(self):
        inp = ['19viewid { Q277   A}']
        ref = {
                'view_num': 19,
                'flipped': False,
                'resn': 'Q',
                'resi': '277',
                'alt': '',
                'chain': 'A',
                }
        res = kin.process_viewid(inp, self.base_context)
        self.assertEqual(res, ref)

    def test_with_flipped_asterisk(self):
        inp = ['2viewid {*Q32   A}']
        res = kin.process_viewid(inp, self.base_context)
        self.assertTrue(res['flipped'])

    def test_with_insertion_code(self):
        inp = ['2viewid { Q32A  A}']
        res = kin.process_viewid(inp, self.base_context)
        self.assertEqual(res['resi'], '32A')

    @mock.patch('pymolprobity.kinemage.logger')
    def test_with_bad_format(self, mock_logger):
        inp = ['viewid bad bad bad']
        res = kin.process_viewid(inp, self.base_context)
        self.assertTrue(mock_logger.warning.called)
        self.assertIsNone(res)


@mock.patch('pymolprobity.kinemage.single_line_keyword_check')
class ProcessMasterTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_master(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_well_formed_master(self, mock_check):
        inp = ['master {something}']
        res = kin.process_master(inp, self.base_context)
        self.assertEqual(res, 'something')


@mock.patch('pymolprobity.kinemage.single_line_keyword_check')
class ProcessPointmasterTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_pointmaster(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_well_formed_pointmaster(self, mock_check):
        inp = ["pointmaster 'a' {something}"]
        res = kin.process_pointmaster(inp, self.base_context)
        ref = {'code': 'a', 'label': 'something', 'enable': 1}
        self.assertEqual(res, ref)

    def test_with_on_statement(self, mock_check):
        inp = ["pointmaster 'a' {something} on"]
        res = kin.process_pointmaster(inp, self.base_context)
        ref = {'code': 'a', 'label': 'something', 'enable': 1}
        self.assertEqual(res, ref)

    def test_with_off_statement(self, mock_check):
        inp = ["pointmaster 'a' {something} off"]
        res = kin.process_pointmaster(inp, self.base_context)
        ref = {'code': 'a', 'label': 'something', 'enable': 0}
        self.assertEqual(res, ref)


@mock.patch('pymolprobity.kinemage.single_line_keyword_check')
class ProcessKinemageKeywordTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_kinemage_keyword(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_well_formed_kinemage(self, mock_check):
        inp = ['kinemage 1']
        res = kin.process_kinemage_keyword(inp, self.base_context)
        self.assertEqual(res, '1')


@mock.patch('pymolprobity.kinemage.single_line_keyword_check')
class ProcessGroupTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_group(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_dominant_group(self, mock_check):
        inp = ['group {something} dominant']
        res = kin.process_group(inp, self.base_context)
        ref = ['something', 'dominant']
        self.assertEqual(res, ref)

    def test_with_animate_group(self, mock_check):
        inp = ['group {something} animate']
        res = kin.process_group(inp, self.base_context)
        ref = ['something', 'animate']
        self.assertEqual(res, ref)


@mock.patch('pymolprobity.kinemage.single_line_keyword_check')
class ProcessSubgroupTests(unittest.TestCase):
    def setUp(self):
        self.base_context = {
                'kinemage': None,
                'group': None,
                'subgroup': None,
                'animate': 0,
            }

    def test_calls_single_line_keyword_check(self, mock_check):
        inp = ['blah']
        res = kin.process_subgroup(inp, self.base_context)
        self.assertTrue(mock_check.called)

    def test_with_dominant_subgroup(self, mock_check):
        inp = ['subgroup {something} dominant']
        res = kin.process_subgroup(inp, self.base_context)
        ref = [None, 'something', 'dominant', None]
        self.assertEqual(res, ref)

    def test_with_dominant_before_name_subgroup(self, mock_check):
        inp = ['subgroup dominant {something}']
        res = kin.process_subgroup(inp, self.base_context)
        ref = ['dominant', 'something', None, None]
        self.assertEqual(res, ref)

    def test_with_nobutton_dominant_subgroup(self, mock_check):
        inp = ['subgroup {something} nobutton dominant']
        res = kin.process_subgroup(inp, self.base_context)
        ref = [None, 'something', 'nobutton', 'dominant']
        self.assertEqual(res, ref)



