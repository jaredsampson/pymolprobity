'''Basic tests for PyMOLProbity plugin.'''

import mock
import unittest

# PyMOL API setup
import __main__
__main__.pymol_argv = ['pymol','-qkc']
import pymol
from pymol import cmd

from .context import pymolprobity
import pymolprobity.main as mp
import pymolprobity.flips as flp


###############################################################################
#
#  GENERAL FUNCTION TESTS
#
###############################################################################

class GetObjectTests(unittest.TestCase):
    def setUp(self):
        '''Create an MPObject and register it with the module-level objects dict.'''
        self.obj = mp.MPObject('test')
        mp.objects['test'] = self.obj

    def tearDown(self):
        del(self.obj)
        mp.objects.clear()

    def test_string_existing_object(self):
        '''Retrieve the MPObject by its string name.'''
        o = mp.get_object('test')
        assert o == self.obj
        assert o.name == 'test'

    def test_string_nonexistent_object(self):
        '''A string not in the objects dict returns None.'''
        with self.assertRaises(AttributeError):
            mp.get_objects('blah')

    def test_mpobject_in_objects_dict(self):
        '''Retrieve the MPObject itself.'''
        assert mp.get_object(self.obj) == self.obj

    def test_mpobject_not_in_objects_dict(self):
        '''An MPObject not in the module level objects dict also returns None.'''
        other_obj = mp.MPObject('other')
        with self.assertRaises(ValueError):
            mp.get_object(other_obj)

    def test_other_input_type(self):
        '''Non-string and non-MPObject input returns None.'''
        with self.assertRaises(ValueError):
            mp.get_object(1)


@mock.patch('pymolprobity.main.MPObject')
@mock.patch('pymolprobity.main.get_object')
class GetOrCreateObjectTests(unittest.TestCase):
    def setUp(self):
        # should probably mock out main.objects dict instead
        mp.objects['test'] = mp.MPObject('test')

    def tearDown(self):
        mp.objects.clear()

    def test_get_existing_object(self, mock_get, mock_MPObj):
        res = mp.get_or_create_object('test')
        ref = mock_get.return_value
        self.assertEqual(res, ref)
        mock_MPObj.assert_not_called()

    def test_with_nonexisting_object(self, mock_get, mock_MPObj):
        res = mp.get_or_create_object('blah')
        ref = mock_get.return_value
        self.assertEqual(res, ref)
        mock_get.assert_called_once_with('blah')

    def test_non_string_input(self, mock_get, mock_MPObj):
        with self.assertRaises(TypeError):
            mp.get_or_create_object(2)


class RunCommandTests(unittest.TestCase):
    def test_with_args_only(self):
        output = mp.run_command(['printf', 'blah'])
        assert output == "blah"

    def test_with_args_and_input_str(self):
        output = mp.run_command(['cat', '-'], 'blah')
        assert output == "blah"

    @mock.patch('pymolprobity.main.logger')
    def test_nonzero_return_code(self, mock_logger):
        # TODO mock subprocess
        mp.run_command(['false'])
        self.assertTrue(mock_logger.warning.called)

    @mock.patch('pymolprobity.main.logger')
    def test_with_missing_executable(self, mock_logger):
        mp.run_command(['not_a_command'])
        self.assertTrue(mock_logger.error.called)


@mock.patch('pymolprobity.main.cmd.save')
@mock.patch('pymolprobity.main.tempfile')
class SaveToTempfileTests(unittest.TestCase):
    def test_basic(self, mock_tempfile, mock_save):
        inp = 'mystr'
        fn = 'somefile.pdb'
        tf = mock_tempfile.NamedTemporaryFile.return_value
        tf.name = 'somefile.pdb'
        res = mp.save_to_tempfile(inp)
        tf.write.assert_called_once_with(inp)
        self.assertEqual(res, fn)



###############################################################################
#
#  MPOBJECT TESTS
#
###############################################################################

class MPObjectTests(unittest.TestCase):
    def setUp(self):
        self.o = mp.MPObject('test')

    @mock.patch('pymolprobity.main.cmd')
    def test_get_pdbstr_with_defaults(self, mock_cmd):
        mock_cmd.get_names.return_value = ['mp_test.test_reduce']
        mock_cmd.get_pdbstr.return_value = 'pdbstr'
        res = self.o.get_pdbstr()
        self.assertEqual(res, 'pdbstr')
        mock_cmd.get_pdbstr.assert_called_once_with('mp_test.test_reduce')

    @mock.patch('pymolprobity.main.cmd')
    def test_get_pdbstr_with_explicit_reduce(self, mock_cmd):
        mock_cmd.get_names.return_value = ['mp_test.test_reduce']
        mock_cmd.get_pdbstr.return_value = 'pdbstr'
        res = self.o.get_pdbstr('reduce')
        self.assertEqual(res, 'pdbstr')
        mock_cmd.get_pdbstr.assert_called_once_with('mp_test.test_reduce')

    @mock.patch('pymolprobity.main.cmd')
    def test_get_pdbstr_with_missing_object(self, mock_cmd):
        mock_cmd.get_names.return_value = ['not here']
        res = self.o.get_pdbstr()
        mock_cmd.get_pdbstr.assert_not_called()
        self.assertEqual(res, None)

    def test_get_kin_cgo_group_with_valid_kin(self):
        res = self.o.get_kin_cgo_group('flipkinNQ')
        ref = 'mp_test.flipkinNQ'
        self.assertEqual(res, ref)

    @mock.patch('pymolprobity.main.logger')
    def test_get_kin_cgo_group_with_invalid_kin(self, mock_logger):
        res = self.o.get_kin_cgo_group('reduce')
        ref = None
        self.assertEqual(res, ref)
        self.assertTrue(mock_logger.warning.called)

    # @mock.patch('pymolprobity.main.cmd')
    # def test_get_kin_group_cgo_group_with_valid_kin(self, mock_cmd):
    #     ref = 'mp_test.flipkinNQ.reduce'
    #     mock_cmd.get_names.return_value = [ref]
    #     res = self.o.get_kin_group_cgo_group('flipkinNQ', 'reduce')
    #     self.assertEqual(res, ref)

    # @mock.patch('pymolprobity.main.logger')
    # @mock.patch('pymolprobity.main.cmd')
    # def test_get_kin_group_cgo_group_with_invalid_kin(self, mock_cmd, mock_logger):
    #     mock_cmd.get_names.return_value = []
    #     res = self.o.get_kin_group_cgo_group('flipkinNQ', 'blah')
    #     ref = None
    #     self.assertEqual(res, ref)
    #     self.assertTrue(mock_logger.warning.called)

    @mock.patch('pymolprobity.main.MPObject.solo_pdb')
    @mock.patch('pymolprobity.main.MPObject.solo_kin')
    @mock.patch('pymolprobity.main.MPObject.get_kin_cgo_group')
    @mock.patch('pymolprobity.main.kinemage.Kinemage', autospec=True)
    @mock.patch('pymolprobity.main.MPObject.disable_kin')
    def test_draw_with_no_args(self, mock_disable_kin, mock_kin, mock_get_grp,
            mock_solo_kin, mock_solo_pdb):
        k = 'probe'
        self.o.kin[k] = mock_kin.return_value
        self.o.draw(k)
        mock_disable_kin.assert_called_once_with()
        mock_get_grp.assert_called_once_with(k)
        grp = mock_get_grp.return_value
        self.o.kin[k].draw.assert_called_with(grp, dot_mode=0)
        mock_solo_pdb.assert_called_once_with(k)
        mock_solo_kin.assert_called_once_with(k)

    @mock.patch('pymolprobity.main.MPObject.animate')
    @mock.patch('pymolprobity.main.MPObject.solo_pdb')
    @mock.patch('pymolprobity.main.MPObject.solo_kin')
    @mock.patch('pymolprobity.main.MPObject.get_kin_cgo_group')
    @mock.patch('pymolprobity.main.kinemage.Kinemage', autospec=True)
    @mock.patch('pymolprobity.main.MPObject.disable_kin')
    def test_draw_with_flipkin_arg(self, mock_disable_kin, mock_kin, mock_get_grp,
            mock_solo_kin, mock_solo_pdb, mock_animate):
        k = 'flipkinNQ'
        self.o.kin[k] = mock_kin.return_value
        self.o.draw(k)
        grp = mock_get_grp.return_value
        self.assertTrue(self.o.kin[k].draw.called)
        self.assertTrue(mock_animate.called)

    @mock.patch('pymolprobity.main.cmd')
    def test_disable_pdb_with_default_all_input(self, mock_cmd):
        for p in self.o.pdb.keys():
            self.o.pdb[p] = '{}_pdb_obj'.format(p)
        self.o.disable_pdb()  # default is pdb='all'
        mock_cmd.disable.assert_has_calls(
                [mock.call('reduce_pdb_obj'),
                 mock.call('flipkinNQ_pdb_obj'),
                 mock.call('flipkinH_pdb_obj'),
                 mock.call('probe_pdb_obj')], any_order=True)

    @mock.patch('pymolprobity.main.cmd')
    def test_disable_pdb_with_valid_single_input(self, mock_cmd):
        inp = 'reduce'
        self.o.pdb[inp] = 'pdb_obj'
        self.o.disable_pdb(inp)
        mock_cmd.disable.assert_called_once_with('pdb_obj')

    @mock.patch('pymolprobity.main.cmd')
    def test_disable_pdb_with_invalid_single_input(self, mock_cmd):
        inp = 'blah'
        with self.assertRaises(KeyError):
            self.o.disable_pdb(inp)

    @mock.patch('pymolprobity.main.cmd')
    def test_enable_pdb(self, mock_cmd):
        inp = 'reduce'
        self.o.pdb[inp] = 'obj'
        self.o.enable_pdb(inp)
        mock_cmd.enable.assert_called_once_with('obj')

    @mock.patch('pymolprobity.main.MPObject.enable_pdb')
    @mock.patch('pymolprobity.main.MPObject.disable_pdb')
    def test_solo_pdb(self, mock_disable, mock_enable):
        inp = 'pdb'
        self.o.solo_pdb(inp)
        mock_disable.assert_called_once_with('all')
        mock_enable.assert_called_once_with(inp)

    @mock.patch('pymolprobity.main.cmd')
    @mock.patch('pymolprobity.main.MPObject.get_kin_cgo_group', autospec=True)
    def test_disable_kin_with_default(self, mock_get_grp, mock_cmd):
        grp = mock_get_grp.return_value
        self.o.disable_kin()
        mock_get_grp.assert_has_calls(
                [mock.call(self.o, 'flipkinNQ'),
                 mock.call(self.o, 'flipkinH'),
                 mock.call(self.o, 'probe')], any_order=True)
        mock_cmd.disable.assert_has_calls(
                [mock.call(grp),
                 mock.call(grp),
                 mock.call(grp)])

    @mock.patch('pymolprobity.main.cmd')
    @mock.patch('pymolprobity.main.MPObject.get_kin_cgo_group', autospec=True)
    def test_disable_kin_with_kin_key(self, mock_get_grp, mock_cmd):
        inp = 'flipkinNQ'
        grp = mock_get_grp.return_value
        self.o.disable_kin(inp)
        mock_get_grp.assert_called_once_with(self.o, inp)
        mock_cmd.disable.assert_called_once_with(grp)

    @mock.patch('pymolprobity.main.cmd')
    @mock.patch('pymolprobity.main.MPObject.get_kin_cgo_group', autospec=True)
    def test_enable_kin_with_kin_key(self, mock_get_grp, mock_cmd):
        inp = 'flipkinNQ'
        grp = mock_get_grp.return_value
        self.o.enable_kin(inp)
        mock_get_grp.assert_called_once_with(self.o, inp)
        mock_cmd.enable.assert_called_once_with(grp)

    @mock.patch('pymolprobity.main.MPObject.enable_kin', autospec=True)
    @mock.patch('pymolprobity.main.MPObject.disable_kin', autospec=True)
    def test_solo_kin_with_kin_key(self, mock_disable, mock_enable):
        inp = 'flipkinNQ'
        self.o.solo_kin(inp)
        mock_disable.assert_called_once_with(self.o, 'all')
        mock_enable.assert_called_once_with(self.o, inp)

    @mock.patch('pymolprobity.main.cmd')
    def test_enable_flipkin_group_with_default_reduce(self, mock_cmd):
        self.o.enable_flipkin_group()
        ref = 'mp_test.*.reduce'
        mock_cmd.enable.assert_called_once_with(ref)













    # @mock.patch('pymolprobity.main.MPObject.get_vectors_cgo')
    # @mock.patch('pymolprobity.main.MPObject.get_dots_cgo')
    # @mock.patch('pymolprobity.main.cmd')
    # def test_draw_workflow_with_defaults(self, mock_cmd, mock_dots_cgo,
    #         mock_vectors_cgo):
    #     mock_cmd.get_view.return_value = 'view'
    #     dots_cgo = mock_dots_cgo.return_value
    #     vectors_cgo = mock_vectors_cgo.return_value

    #     self.o.draw()

    #     mock_cmd.get_view.assert_called_once_with()
    #     mock_cmd.set.assert_has_calls(
    #             [mock.call('cgo_use_shader',     0),
    #              mock.call('cgo_sphere_quality', 0)])
    #     mock_dots_cgo.assert_called_once_with(dot_mode=0)
    #     mock_vectors_cgo.assert_called_once_with()
    #     mock_cmd.load_cgo.assert_has_calls(
    #             [mock.call(dots_cgo,    'test_dots'   ),
    #              mock.call(vectors_cgo, 'test_clashes')])
    #     mock_cmd.group.assert_has_calls(
    #             [mock.call('mp_test', 'test_dots'),
    #              mock.call('mp_test', 'test_clashes')])
    #     mock_cmd.set_view.assert_called_once_with('view')

    # @mock.patch('pymolprobity.main.MPObject.get_vectors_cgo')
    # @mock.patch('pymolprobity.main.MPObject.get_dots_cgo')
    # @mock.patch('pymolprobity.main.cmd')
    # def test_draw_with_dot_mode_1(self, mock_cmd, mock_dots_cgo,
    #         mock_vectors_cgo):
    #     dots_cgo = mock_dots_cgo.return_value
    #     self.o.draw(dot_mode=1)
    #     mock_dots_cgo.assert_called_once_with(dot_mode=1)

    # @mock.patch('pymolprobity.main.points.Dot')
    # def test_get_dots_cgo(self, mock_dot):
    #     o = mp.MPObject('test')
    #     o.dots = [mock_dot, mock_dot, mock_dot]
    #     mock_dot.get_cgo.return_value = ['dots_cgo']
    #     res = o.get_dots_cgo()
    #     ref = ['dots_cgo', 'dots_cgo', 'dots_cgo']
    #     self.assertEqual(res, ref)

    # @mock.patch('pymolprobity.main.points.Vector')
    # def test_get_vectors_cgo(self, mock_vec):
    #     o = mp.MPObject('test')
    #     o.vectors = [mock_vec, mock_vec, mock_vec]
    #     mock_vec.get_cgo.return_value = ['vectors_cgo']
    #     res = o.get_vectors_cgo()
    #     ref = ['vectors_cgo', 'vectors_cgo', 'vectors_cgo']
    #     self.assertEqual(res, ref)

    # def test_get_flip_matching_atom(self):
    #     f0 = flp.Flip()
    #     f0.chain, f0.resi, f0.resn, f0.alt = ['A', '1', 'ASN', '']
    #     f1 = flp.Flip()
    #     f1.chain, f1.resi, f1.resn, f1.alt = ['A', '100', 'GLN', '']

    #     o = mp.MPObject('test')
    #     o.flips = [f0, f1]

    #     atom = {'chain': 'A', 'resi': '100', 'resn': 'GLN', 'alt': ''}
    #     res = o.get_flip_matching_atom(atom)
    #     self.assertEqual(res, f1)

    # @mock.patch('pymolprobity.main.logger')
    # def test_get_flip_matching_atom_with_no_match(self, mock_logger):
    #     f0 = flp.Flip()
    #     f0.chain, f0.resi, f0.resn, f0.alt = ['A', '1', 'ASN', '']
    #     f1 = flp.Flip()
    #     f1.chain, f1.resi, f1.resn, f1.alt = ['A', '100', 'GLN', '']

    #     o = mp.MPObject('test')
    #     o.flips = [f0, f1]

    #     atom = {'chain': 'A', 'resi': '2', 'resn': 'GLN', 'alt': ''}
    #     res = o.get_flip_matching_atom(atom)
    #     self.assertEqual(res, None)

    # @mock.patch('pymolprobity.main.logger')
    # def test_get_flip_matching_atom_with_multiple_matches(self, mock_logger):
    #     f0 = flp.Flip()
    #     f0.chain, f0.resi, f0.resn, f0.alt = ['A', '100', 'ASN', '']

    #     o = mp.MPObject('test')
    #     o.flips = [f0, f0]

    #     atom = {'chain': 'A', 'resi': '100', 'resn': 'ASN', 'alt': ''}
    #     res = o.get_flip_matching_atom(atom)
    #     self.assertEqual(res, None)







###############################################################################
#
#  REDUCE TESTS
#
###############################################################################

class GetReduceArgsTests(unittest.TestCase):
    def test_default_values(self):
        '''Default arguments'''
        ref = ['reduce', '-Quiet', '-FLIP', '-']
        args = mp.get_reduce_args()
        assert args == ref

    def test_no_h(self):
        ref = ['reduce', '-Quiet', '-FLIP', '-Trim', '-']
        args = mp.get_reduce_args(h=0)
        assert args == ref

    def test_no_flips(self):
        ref = ['reduce', '-Quiet', '-NOFLIP', '-']
        args = mp.get_reduce_args(flip=0)
        assert args == ref

    def test_no_quiet(self):
        '''Without -Quiet flag.'''
        ref = ['reduce', '-FLIP', '-']
        args = mp.get_reduce_args(quiet=0)
        assert args == ref

    def test_with_addflags_string(self):
        ref = ['reduce', '-Quiet', '-FLIP', 'added', '-']
        args = mp.get_reduce_args(addflags='added')
        assert args == ref

    def test_with_h_nonint(self):
        with self.assertRaises(TypeError):
            mp.get_reduce_args(h='str')

    def test_with_h_nonzero_int(self):
        ref = ['reduce', '-Quiet', '-FLIP', '-']
        args = mp.get_reduce_args(h=2)

    def test_with_flip_nonint(self):
        with self.assertRaises(TypeError):
            mp.get_reduce_args(flip='str')

    def test_with_flip_nonzero_int(self):
        ref = ['reduce', '-Quiet', '-FLIP', '-']
        args = mp.get_reduce_args(flip=2)

    def test_with_quiet_nonint(self):
        with self.assertRaises(TypeError):
            mp.get_reduce_args(quiet='str')

    def test_with_quiet_nonzero_int(self):
        ref = ['reduce', '-Quiet', '-FLIP', '-']
        args = mp.get_reduce_args(quiet=2)

    def test_with_addflags_nonstring(self):
        with self.assertRaises(TypeError):
            mp.get_reduce_args(addflags=1)


@mock.patch('pymolprobity.main.run_command')
@mock.patch('pymolprobity.main.get_reduce_args')
class GenerateReduceResultTests(unittest.TestCase):

    def test_with_flip_type_0(self, mock_args, mock_run):
        mp.generate_reduce_output('pdbstr', 0)
        mock_args.assert_called_once_with(flip=0)

    def test_with_flip_type_1(self, mock_args, mock_run):
        mp.generate_reduce_output('pdbstr', 1)
        mock_args.assert_called_once_with()

    def test_with_flip_type_2(self, mock_args, mock_run):
        mp.generate_reduce_output('pdbstr', 2)
        mock_args.assert_called_once_with(addflags='-NOBUILD0')

    def test_run_reduce_call(self, mock_args, mock_run):
        mock_args.return_value = ['args']
        mp.generate_reduce_output('pdbstr', 1)
        mock_run.assert_called_once_with(['args'], 'pdbstr')


@mock.patch('pymolprobity.main.logger')
@mock.patch('pymolprobity.main.process_reduce_output')
@mock.patch('pymolprobity.main.generate_reduce_output')
@mock.patch('pymolprobity.main.cmd')
class ReduceObjectTests(unittest.TestCase):
    def test_default_with_flips(self, mock_cmd, mock_gen, mock_proc,
            mock_logger):
        mock_cmd.get_pdbstr.return_value = 'pdbstr'
        mock_gen.return_value = 'reduced_pdbstr'
        mock_proc.return_value = ['flip1', 'flip2']
        mp.reduce_object('obj')
        mock_cmd.get_pdbstr.assert_called_once_with('obj')
        mock_gen.assert_called_once_with('pdbstr', flip_type=1)
        mock_proc.assert_called_once_with('reduced_pdbstr')

    def test_without_flips(self, mock_cmd, mock_gen, mock_proc, mock_logger):
        mock_cmd.get_pdbstr.return_value = 'pdbstr'
        mp.reduce_object('obj', flip=0)
        mock_cmd.get_pdbstr.assert_called_once_with('obj')
        mock_gen.assert_called_once_with('pdbstr', flip_type=0)

    def test_no_reduce_output(self, mock_cmd, mock_gen, mock_proc,
            mock_logger):
        mock_gen.return_value = None
        mp.reduce_object('obj')
        self.assertTrue(mock_logger.error.called)
        self.assertFalse(mock_proc.called)


@mock.patch('pymolprobity.flips.parse_flips')
class ProcessReduceOutputTests(unittest.TestCase):
    def test_with_normal_string_input(self, mock_parse_flips):
        val = 'USER  MOD blah blah\nHEADER blah\nUSER  MOD blah blah blah'
        user_mod_ref = ['blah blah', 'blah blah blah']
        mp.process_reduce_output(val)
        mock_parse_flips.assert_called_with(user_mod_ref)

    def test_with_nonstring_input(self, mock_parse):
        '''Raise a TypeError with non-string input.'''
        bad_val = 12345
        with self.assertRaises(TypeError):
            mp.process_reduce_output(bad_val)

    def test_with_malformed_string_input(self, mock_parse):
        '''Raise a ValueError if the input string doesn't look like normal
        Reduce output.'''
        bad_val = "This doesn't start with USER  MOD."
        with self.assertRaises(ValueError):
            mp.process_reduce_output(bad_val)

    def test_calls_parse_flips(self, mock_parse):
        val = "USER  MOD blah blah\nHEADER blah\nUSER  MOD blah blah blah"
        user_mod_ref = ["blah blah", "blah blah blah"]
        mp.process_reduce_output(val)
        mock_parse.assert_called_once_with(user_mod_ref)


###############################################################################
#
#  FLIPKIN TESTS
#
###############################################################################

@mock.patch('pymolprobity.main.logger')
@mock.patch('pymolprobity.main.os.unlink')
@mock.patch('pymolprobity.main.process_flipkin_output')
@mock.patch('pymolprobity.main.generate_flipkin_output')
@mock.patch('pymolprobity.main.save_to_tempfile')
@mock.patch('pymolprobity.main.get_object')
class FlipkinObjectTests(unittest.TestCase):
    def test_workflow_for_nq_defaults(self, mock_get_obj, mock_tf, mock_gen,
            mock_proc, mock_unlink, mock_logger):  #, mock_apply):
        o = mock_get_obj.return_value
        o.reduce_output = 'pdbstr'
        mock_tf.return_value = 'tempfile.pdb'
        mock_gen.side_effect = ['flipkin_nq', 'flipkin_h']
        mock_proc.side_effect = ['proc_nq', 'proc_h']

        mp.flipkin_object('test')

        mock_get_obj.assert_called_once_with('test')
        mock_tf.assert_called_once_with('pdbstr')
        mock_gen.assert_has_calls(
                [mock.call('tempfile.pdb'),
                 mock.call('tempfile.pdb', his=True)])
        mock_proc.assert_has_calls(
                [mock.call('flipkin_nq'),
                 mock.call('flipkin_h')])
        mock_unlink.assert_called_once_with('tempfile.pdb')
        # mock_apply.assert_has_calls(
        #         [mock.call(o, 'proc_nq'),
        #          mock.call(o, 'proc_h')])

    def test_no_flipkinNQ_output(self, mock_get_obj, mock_tf, mock_gen,
            mock_proc, mock_unlink, mock_logger):
        mock_gen.side_effect = (None, 'flipkin_h')
        mp.flipkin_object('test')
        self.assertTrue(mock_logger.error.called)
        self.assertFalse(mock_proc.called)

    def test_no_flipkinH_output(self, mock_get_obj, mock_tf, mock_gen,
            mock_proc, mock_unlink, mock_logger):
        mock_gen.side_effect = ('flipkin_nq', None)
        mp.flipkin_object('test')
        self.assertTrue(mock_logger.error.called)
        self.assertFalse(mock_proc.called)




@mock.patch('pymolprobity.main.run_command')
class GenerateFlipkinOutputTests(unittest.TestCase):
    def test_workflow_with_default_nq(self, mock_run):
        mock_run.return_value = 'output'
        res = mp.generate_flipkin_output('filename')
        mock_run.assert_called_once_with(['flipkin', 'filename'])
        ref = 'output'
        self.assertEqual(res, ref)

    def test_workflow_with_his(self, mock_run):
        mp.generate_flipkin_output('filename', his=True)
        mock_run.assert_called_once_with(['flipkin', '-h', 'filename'])


@mock.patch('pymolprobity.main.kinemage.process_kinemage')
class ProcessFlipkinOutputTests(unittest.TestCase):
    def test_workflow(self, mock_proc_kin):
        kin = mock_proc_kin.return_value
        res = mp.process_flipkin_output('kinstr')
        mock_proc_kin.assert_called_once_with('kinstr')


# @mock.patch('pymolprobity.points.Vector')
# @mock.patch('pymolprobity.flips.Flip')
# @mock.patch('pymolprobity.main.kinemage.Kinemage')
# @mock.patch('pymolprobity.main.MPObject')
# class ApplyFlipkinToMPObject(unittest.TestCase):
#     def test_with_animated_coordinate_vectorlist_input(self, mock_Obj, mock_Kin,
#             mock_Flip, mock_Vector):
#         mpobj = mock_Obj('myobj')
#         f = mock_Flip()
#         f.flipped = False
#         mpobj.get_flip_matching_atom.return_value = f
#         kin = mock_Kin()
#         vec = mock_Vector()
#         vec.atom = [{'name': 'ca', 'coords': [1, 2, 3], 'resn': 'GLN'},
#                     {'name': 'cb', 'coords': [2, 3, 4], 'resn': 'ASN'}]
#         kin.keywords = {
#                 '1': {'keyword': 'group', 'data': ['reduce', 'animate']},
#                 '2': {'keyword': 'vectorlist', 'data': [vec]},
#                 }
#         mp.apply_flipkin_to_mpobject(mpobj, kin)
#         mpobj.get_flip_matching_atom.assert_called_once_with(vec.atom[1])

#     def test_with_coordinates_vectorlist_input(self, mock_Obj, mock_Kin,
#             mock_Flip, mock_Vector):
#         '''should ignore vectorlist within main coordinates group'''
#         mpobj = mock_Obj('myobj')
#         kin = mock_Kin()
#         vec = mock_Vector()
#         kin.keywords = {
#                 '1': {'keyword': 'group', 'data': ['objname', 'dominant']},
#                 '2': {'keyword': 'vectorlist', 'data': [vec]},
#                 }
#         mp.apply_flipkin_to_mpobject(mpobj, kin)
#         self.assertFalse(mpobj.get_flip_matching_atom.called)

#     @mock.patch('pymolprobity.main.logger')
#     def test_with_kinemage_keyword_input(self, mock_logger, mock_Obj,
#             mock_Kin, mock_Flip, mock_Vector):
#         mpobj = mock_Obj('myobj')
#         kin = mock_Kin()
#         vec = mock_Vector()
#         kin.keywords = {
#                 '1': {'keyword': 'kinemage', 'data': '1'}}
#         mp.apply_flipkin_to_mpobject(mpobj, kin)
#         self.assertTrue(mock_logger.debug.called)

#     @mock.patch('pymolprobity.main.logger')
#     def test_with_subgroup_input(self, mock_logger, mock_Obj,
#             mock_Kin, mock_Flip, mock_Vector):
#         mpobj = mock_Obj('myobj')
#         kin = mock_Kin()
#         vec = mock_Vector()
#         kin.keywords = {
#                 '1': {'keyword': 'subgroup', 'data': ['foo', 'bar']}}
#         mp.apply_flipkin_to_mpobject(mpobj, kin)
#         mock_logger.debug.assert_called_with('entering subgroup: bar')






###############################################################################
#
#  PROBE TESTS
#
###############################################################################

class GetProbeArgsTests(unittest.TestCase):
    def test(self):
        fn = 'temp.pdb'
        ref = ['probe', '-Quiet', '-Self', 'ALL', 'temp.pdb']
        res = mp.get_probe_args(fn)
        self.assertEqual(res, ref)


@mock.patch('pymolprobity.main.run_command')
@mock.patch('pymolprobity.main.get_probe_args')
@mock.patch('pymolprobity.main.os')
@mock.patch('pymolprobity.main.save_to_tempfile')
class GenerateProbeOutputTests(unittest.TestCase):
    def test_workflow(self, mock_tf, mock_os, mock_args,
            mock_run):
        fn = 'somefile.pdb'
        pdb = 'pdbstr'
        mock_tf.return_value = fn
        mock_os.path.isfile.side_effect = [True, False]
        mock_args.return_value = [1,2,3]
        mock_run.return_value = 'output'

        res = mp.generate_probe_output(pdb)
        self.assertEqual(res, 'output')

        mock_tf.assert_called_once_with(pdb)
        self.assertTrue(mock_args.called)
        mock_run.assert_called_once_with([1, 2, 3])
        mock_os.unlink.assert_called_once_with(fn)
        call = mock.call(fn)
        mock_os.path.isfile.assert_has_calls([call, call])


@mock.patch('pymolprobity.main.logger')
@mock.patch('pymolprobity.main.cmd')
@mock.patch('pymolprobity.main.process_probe_output')
@mock.patch('pymolprobity.main.generate_probe_output')
@mock.patch('pymolprobity.main.get_object')
class ProbeObjectTests(unittest.TestCase):
    def test_workflow(self, mock_get_obj, mock_gen, mock_proc,
            mock_cmd, mock_logger):
        obj = 'obj'
        o = mock_get_obj.return_value
        # o.pdb['userflips'] == 'userflips_obj'
        # o.pdb['reduce'] == 'reduce_obj'
        o.get_pdbstr.return_value = 'pdbstr'
        mock_gen.return_value = 'output'
        mock_proc.return_value = ('dots_list', 'clashes_list')

        mp.probe_object(obj)

        mock_get_obj.assert_called_once_with(obj)
        mock_gen.assert_called_once_with('pdbstr')
        self.assertTrue(o.draw.called)

    def test_no_probe_output(self, mock_get_obj, mock_gen, mock_proc,
            mock_cmd, mock_logger):
        '''Fail gracefully if `probe` call returns no output.'''
        mock_gen.return_value = None
        mp.probe_object('obj')
        self.assertTrue(mock_logger.error.called)
        self.assertFalse(mock_proc.called)




class ProcessProbeOutputTests(unittest.TestCase):
    @mock.patch('pymolprobity.main.kinemage.process_kinemage')
    def test_workflow(self, mock_proc_kin):
        kin = mock_proc_kin.return_value
        inp = '@some\n@kinemage\n@string'
        res = mp.process_probe_output(inp)
        mock_proc_kin.assert_called_once_with(inp)
        self.assertEqual(res, kin)





if __name__ == '__main__':
    unittest.main()