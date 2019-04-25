'''Main functions for PyMOLProbity plugin.'''

from __future__ import absolute_import
from __future__ import print_function

import logging
import os
import re
import subprocess
import tempfile

from pymol import cmd  #cgo, cmd, plugins
from pymol import CmdException

# import colors
# from commands import command_line_output, save_to_tempfile
from . import flips
from . import kinemage
from . import points
# import settings
# import utils


# Set up logger
logger = logging.getLogger(__name__)


# # Load saved settings from ~/.pymolpluginsrc if available
# settings.load_settings()



###############################################################################
#
# EXCEPTIONS
#
###############################################################################

class MPException(CmdException):
    pass


###############################################################################
#
#  MPOBJECT
#
###############################################################################

class MPObject(object):
    '''Container for PyMOLProbity data associated with a loaded PyMOL object.'''

    def get_kin_cgo_group(self, kin):
        '''Return the name of the PyMOL group containing the kinemage CGO objects.

        PARAMETERS

            kin         (str) Name of the kinemage. Possible values are
                        'flipkinNQ', 'flipkinH', or 'probe'.

        '''
        if kin in self.kin.keys():
            return '{}.{}'.format(self.mp_group, kin)
        else:
            msg = 'Unexpected kinemage name specified: {}'.format(kin)
            logger.warning(msg)
            return None

    def draw(self, kin, dot_mode=0):
        '''Draw kinemage CGO objects for this MPObject.

        PARAMETERS

            kin         (str) The kinemage to draw. Possible values are
                        'reduce',  'flipkinNQ', 'flipkinH', or 'probe'.

            dot_mode    (int) Set the style of the kinemage dots: 0 = spheres,
                        1 = quads.  (Default: 0)

        '''
        # Disable (toggle off) all kinemages from this object
        self.disable_kin()

        # Draw the Kinemage CGO into the proper PyMOL object group
        source_kin = self.kin[kin]
        kin_grp = self.get_kin_cgo_group(kin)
        source_kin.draw(kin_grp, dot_mode=dot_mode)

        # Show only the current kinemage and its corresponding molecule.
        self.solo_pdb(kin)
        self.solo_kin(kin)

        if kin.startswith('flipkin'):
            self.animate()

    def enable_pdb(self, pdb):
        '''Enable an MPObject structure object.

        PARAMETERS

            pdb     (str) The structure to be enabled.  Possible values are
                    'reduce', 'flipkinNQ', 'flipkinH', 'userflips', or 'probe'.

        '''
        pdb_obj = self.pdb[pdb]
        cmd.enable(pdb_obj)

    def disable_pdb(self, pdb='all'):
        '''Disable an MPObject structure object.

        PARAMETERS

            pdb     (str) The structure to be disabled.  Possible values are
                    'reduce', 'flipkinNQ', 'flipkinH', 'userflips', 'probe', or
                    'all'.  (Default: all)

        '''
        if pdb == 'all':
            for p in self.pdb.keys():
                self.disable_pdb(p)
        else:
            pdb_obj = self.pdb[pdb]
            logger.debug('disabling {}'.format(pdb_obj))
            cmd.disable(pdb_obj)

    def solo_pdb(self, pdb):
        '''Enable an MPObject structure object and disable all the others.

        PARAMETERS

            pdb     (str) The structure to be enabled.  Possible values are
                    'reduce', 'flipkinNQ', 'flipkinH', 'userflips', or 'probe'.

        '''
        self.disable_pdb('all')
        self.enable_pdb(pdb)

    def enable_kin(self, kin):
        '''Enable an MPObject kinemage.

        PARAMETERS

            kin     (str) The kinemage to be enabled.  Possible values are
                    'flipkinNQ', 'flipkinH', or 'probe'.

        '''
        kin_group = self.get_kin_cgo_group(kin)
        cmd.enable(kin_group)

    def disable_kin(self, kin='all'):
        '''Disable an MPObject kinemage.

        PARAMETERS

            kin     (str) The kinemage to be disabled.  Possible values are
                    'flipkinNQ', 'flipkinH', 'probe', or 'all'.  (Default: all)

        '''
        if kin == 'all':
            for k in self.kin.keys():
                self.disable_kin(k)
        else:
            kin_grp = self.get_kin_cgo_group(kin)
            logger.debug('disabling'.format(kin_grp))
            cmd.disable(kin_grp)

    def solo_kin(self, kin):
        '''Enable an MPObject kinemage and disable all the others.

        PARAMETERS

            kin     (str) The kinemage to be enabled.  Possible values are
                    'flipkinNQ', 'flipkinH', or 'probe'.

        '''
        self.disable_kin('all')
        self.enable_kin(kin)

    def enable_flipkin_group(self, group='reduce'):
        '''Enable a flipkin kinemage group.

        PARAMETERS

            group   (str) The group to be enabled.  Possible values are
                    'reduce', 'flipNQ', 'flipH', or 'both'.  (Default: both)

        '''
        if group in ['reduce', 'flipNQ', 'flipH']:
            sel = '{}.*.{}'.format(self.mp_group, group)
        else:
            msg = 'not a typical flipkin group name: {}'.format(group)
            logger.warning(msg)
            sel = '{}.*.{}*'.format(self.mp_group, group)
        logger.debug('enabling "{}"'.format(sel))
        cmd.enable(sel)

    def disable_flipkin_group(self, group='all'):
        '''Disable a flipkin kinemage group.

        PARAMETERS

            group   (str) The group to be disabled.  Possible values are
                    'reduce', 'flipNQ', 'flipH', or 'all'.  (Default: both)

        '''
        GROUPS = ['reduce', 'flipNQ', 'flipH']
        if group == 'all':
            # disable each one recursively
            for g in GROUPS:
                self.disable_flipkin_group(g)
        if group not in GROUPS:
            msg = 'not a typical flipkin group name: {}'.format(group)
            logger.warning(msg)
        sel = '{}.*.{}'.format(self.mp_group, group)
        logger.debug('disabling "{}"'.format(sel))
        cmd.disable(sel)

    def solo_flipkin_group(self, grp):
        '''Enable a flipkin kinemage group and disable all the others.

        PARAMETERS

            group   (str) The group to be enabled.  Possible values are
                    'reduce', 'flipNQ', or 'flipH'.

        '''
        self.disable_flipkin_group('all')
        self.enable_flipkin_group(grp)

    def animate(self, no_recurse=False):
        '''Toggle between 'reduce' and 'flip' flipkin kinemage groups.

        The 'flip' group is either 'flipNQ' or 'flipH', depending on which
        kinemage is enabled.

        Note: Using this method is not recommended if you're using the GUI, as
        it doesn't currently affect the state of the flipkin group checkboxes.

        '''
        # TODO integrate GUI animate button with this function.
        # Check which flipkin kinemages are enabled.
        logger.debug('checking which kinemages are enabled')
        enabled_flipkins = []
        for fk in ['flipkinNQ', 'flipkinH']:
            kin_grp = self.get_kin_cgo_group(fk)
            if kin_grp in cmd.get_names(enabled_only=1):
                logger.debug('  {} flipkin is enabled.'.format(fk))
                enabled_flipkins.append(fk)
                continue

        if not enabled_flipkins:
            if no_recurse:
                logger.debug('endless recursion, something went wrong')
                return

            logger.debug('  no flipkin kinemages enabled...enabling flipkinNQ')
            self.solo_kin('flipkinNQ')
            self.animate(True)
            return

        # Regex to match `mp_myobj.*.reduce` (flipkin 'reduce' group)
        reduce_group_regex = re.compile(re.escape(self.mp_group) + r'\.[^\.]+\.reduce$')

        # If at least one flipkin kinemage is enabled AND either the 'reduce'
        # molecule or any 'reduce' kinemage CGO groups are enabled, solo the
        # 'flipX' groups of both flipkin kinemages and the 'flip' molecule of
        # whichever flipkin kinemages are enabled.
        logger.debug('checking if any reduce molecule or kinemage is enabled...')
        reduce_is_enabled = 0
        for name in cmd.get_names(enabled_only=1):
            logger.debug('  checking {} for reduce group match'.format(name))
            if reduce_group_regex.match(name) or name == self.pdb['reduce']:
                logger.debug('  match! {}')
                reduce_is_enabled = 1
                continue
        if reduce_is_enabled:
            logger.debug('reduce was enabled, switching to flips.')
            self.disable_pdb('reduce')
            # Enable the flipkin 'flip' groups
            self.solo_flipkin_group('flip')
            # And the molecules for each enabled flipkin.
            for fk in enabled_flipkins:
                self.enable_pdb(fk)

        # If, on the other hand, no 'reduce' group or molecule is enabled, solo
        # the 'reduce' groups for all kinemages and the 'reduce'
        # coordinates for this MPObject.
        else:
            self.solo_pdb('reduce')
            self.solo_flipkin_group('reduce')

    def get_pdbstr(self, pdb='reduce'):
        '''Return a PDB string for the specified structure object.'''
        # TODO this needs some work
        name = self.pdb[pdb]
        if name in cmd.get_names():
            return cmd.get_pdbstr(name)
        else:
            msg = 'Sorry, no object loaded called {}!'.format(name)
            logger.warning(msg)
            return None

    def __init__(self, name):
        self.name = name
        group = 'mp_{}'.format(name)
        self.mp_group = group

        # Reduce output (need USER MOD headers to pass to flipkin)
        self.reduce_output = None
        self.flips = None

        # Coordinate object names
        self.pdb = {
                'reduce': '{}.{}_reduce'.format(group, name),
                'flipkinNQ': '{}.{}_flipkinNQ'.format(group, name),
                'flipkinH': '{}.{}_flipkinH'.format(group, name),
                'userflips': '{}.{}_userflips'.format(group, name),
                'probe': '{}.{}_probe'.format(group, name),
                }

        # Processed Kinemage instances from Flipkin & Probe
        self.kin = {
                'flipkinNQ': None,
                'flipkinH':  None,
                'probe':  None,
                }

        # Store viewids for persistence in GUI
        self.views = {
                'flipkinNQ': None,
                'flipkinH': None,
                'probe': None,
                }


    def __str__(self):
        return "<MPObject: %s>" % self.name


###############################################################################
#
#  GENERAL FUNCTIONS
#
###############################################################################

def get_object(obj):
    '''Return the matching MPObject instance if it exists.

    PARAMETERS

        obj     (str) The name of a loaded PyMOL structure object.  (Also
                accepts an MPObject instance.)

    '''
    if type(obj) is str:
        try:
            return objects[obj]
        except KeyError:
            msg = "get_object: '{}' not in plugin objects dict.".format(obj)
            logger.error(msg)
            raise MPException(msg)
    elif type(obj) is MPObject:
        try:
            if obj.name not in objects.keys():
                msg = ("get_object: {}'s name attribute not in plugin objects "
                        'dict.').format(obj)
                raise ValueError(msg)
            if objects[obj.name] is obj:
                return obj
            else:
                msg = ('get_object: MPObject {} is not listed in the plugin '
                        'objects dict under {}!').format(obj, obj.name)
                raise ValueError(msg)
        except AttributeError:
            msg = "get_object: '{}' not in plugin objects dict.".format(obj)
            logger.error(msg)
            raise MPException(msg)
    else:
        msg = "get_object: `obj` must be either a string or MPObject instance."
        raise ValueError(msg)


def get_or_create_object(obj):
    '''Return the matching MPObject instance if it exists, or create it.

    PARAMETERS

        obj     (str) Name of a loaded PyMOL structure object.

    '''
    # Check input
    if not type(obj) is str:
        msg = "get_or_create_object: `obj` must be a string"
        raise TypeError(msg)

    # Get the object...
    if obj in objects.keys():
        logger.debug('Using existing MPObject: {}'.format(obj))
        return get_object(obj)
    # Or create it
    else:
        logger.debug('Creating MPObject: {}'.format(obj))
        objects[obj] = MPObject(obj)
        return get_object(obj)


def save_to_tempfile(data_str):
    """Save a selection to a temporary PDB file and return the file name."""
    # text mode (py3 compatibility)
    tf = tempfile.NamedTemporaryFile("w", suffix=".pdb", dir=".", delete=False)
    tf.write(data_str)
    tf.close()
    return tf.name


def run_command(args, input_str=None):
    """Run a command with the given arguments and optional piped STDIN input
    string, and return STDOUT as a string.
    """
    try:
        process = subprocess.Popen(args, stdin=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True, # text mode (py3 compatibility)
                stdout=subprocess.PIPE)
        output, stderr = process.communicate(input_str or "")
    except OSError:
        msg = ("Unable to run the following command:\n\t`{}`\nPlease make "
               "sure {} is installed and can be found on the shell PATH.")
        raise MPException(msg.format(" ".join(args), args[0]))

    logger.debug("===== BEGIN OUTPUT =====\n%s", output)
    logger.debug("===== END OUTPUT =====")

    # Error check
    if process.returncode != 0:
        msg = '{} returned {}'.format(args[0], str(process.returncode))
        logger.warning(msg)

    if stderr:
        logger.warning(stderr)

    return output



###############################################################################
#
#  REDUCE
#
###############################################################################

def get_reduce_args(h=1, flip=1, quiet=1, addflags=None):
    """Return a list of arguments to be used to run Reduce locally.

    Prepare arguments to run Reduce locally using the specified options.  Note
    that the path to the `reduce` executable must be set in the MolProbity
    plugin settings (see [mpset], [mpget]).

    With no keyword arguments, uses Reduce defaults, to perform NQH flips and
    add hydrogens.

    Any additional options passed via `addflags` will take precedence over
    other options specified by the normal keyword arguments.

    USAGE:

        get_reduce_args [h=1, [flip=1, [quiet=1, [addflags=None]]]]

    ARGUMENTS:

        h:
            1 = Add hydrogens.  (default)
            0 = Remove ("trim") hydrogens if present.

        flip:
            If unset, no NQH flips will be performed.  (default=1)

        quiet:
            Suppress reduce's normal console output.  (default=1)

        addflags:
            A space-separated string of additional flags to pass directly to
            Reduce (e.g.  "-DENSITY24 -SHOWSCORE").  (default=None)

    """
    # Check inputs
    msg = "get_reduce_args: `{}` must be of type (or castable as) `{}`"
    try: h = int(h)
    except: raise TypeError(msg.format('h', int))
    try: flip = int(flip)
    except: raise TypeError(msg.format('flip', int))
    try: quiet = int(quiet)
    except: raise TypeError(msg.format('quiet', int))

    # Handle nonbinary values
    if h     != 0:  h     = 1
    if flip  != 0:  flip  = 1
    if quiet != 0:  quiet = 1

    if addflags is not None:
        if type(addflags) is not str:
            msg = ("get_reduce_args: `addflags` should be a space-separated "
                   "string of options to pass to Reduce.")
            raise TypeError(msg)
#     logger.debug('reduce_object inputs are ok.')

    # Begin with the path to `reduce`.
    args = ['reduce']  # [settings.mpgetq('reduce_path')]

    # Handle keyword arguments
    if quiet:
        args.append('-Quiet')

    if flip:
        args.append('-FLIP')
    else:
        args.append('-NOFLIP')

    if h == 0:
        args.append('-Trim')

    # Any user-specified flags are added last.  These will take precedence over
    # any of the previous flags from keyword arguments.
    if addflags:
        args.extend(addflags.split(' '))

    # Pass the pdbstr from STDIN
    args.append('-')

    return args


def generate_reduce_output(pdbstr, flip_type=1):
    '''Generate Reduce output from the given PDB string and flip_type.

    ARGUMENTS
        pdbstr (str)    A loaded PyMOL object (typically an automatically
                        generated temporary one, without hydrogens).

        flip_type (int) Indicates which type of result to make.
                        0 = no flips
                        1 = recommended flips (default)
                        2 = all flips

    '''
    lookup = {
            0: {'flip': 0},               # no flips
            1: {},                        # default flips
            2: {'addflags': '-NOBUILD0'}  # all scorable flips
            }
    args = get_reduce_args(**lookup[flip_type])
    output = run_command(args, pdbstr)
    return output


def process_reduce_output(reduced_pdbstr):
    '''Process Reduce output and return a Reduce result dict.'''
    prefix = 'USER  MOD '
    if not type(reduced_pdbstr) is str:
        msg = 'process_reduce_output: argument `reduced_pdbstr` must be a string.'
        raise TypeError(msg)
    if not reduced_pdbstr.startswith(prefix):
        raise ValueError('process_reduce_output: malformed input string')
    lines = reduced_pdbstr.split('\n')

    # Collect USER MOD records and remove the prefix from each line.
    user_mod = [l.replace(prefix, '') for l in lines if prefix in l]

    flips_list = flips.parse_flips(user_mod)

    return flips_list


def reduce_object(obj, flip=1):
    """Add hydrogens to a copy of a loaded PyMOL object with Reduce.

    TODO: more doc here

    """
    # Run reduce with specified flips
    pdbstr = cmd.get_pdbstr(obj)
    reduced_pdbstr = generate_reduce_output(pdbstr, flip_type=flip)

    # Fail gracefully if no output is generated.
    if not reduced_pdbstr:
        msg = "Failed to generate Reduce output for {}.".format(obj)
        logger.error(msg)
        return

    withflips = " with flips" if flip else ""
    logger.info("Generated Reduce output{} for '{}'.".format(withflips, obj))

    # Process the output string for flips
    flips_list = process_reduce_output(reduced_pdbstr)
    logger.info("Processed Reduce output to extract list of flips.")

    # Store flips list and raw reduced_pdbstr in MPObject
    o = get_or_create_object(obj)
    o.flips = flips_list
    o.reduce_output = reduced_pdbstr

    # Store current group_auto_mode setting
    gam = cmd.get('group_auto_mode')
    cmd.set('group_auto_mode', 2)

    # Load the output PDB into a copy of the original and disable the original.
    name = o.pdb['reduce']
    cmd.create(name, obj)  # duplicate original to preserve representation
    cmd.read_pdbstr(reduced_pdbstr, name, state=1)
    cmd.disable(obj)

    # Restore original group_auto_mode
    cmd.set('group_auto_mode', gam)



###############################################################################
#
#  FLIPKIN
#
###############################################################################

def generate_flipkin_output(filename, his=False):
    flipkin_path = 'flipkin'  # TODO settings
    args = [flipkin_path]
    if his:
        args.append('-h')
    args.append(filename)
    output = run_command(args)
    return output


def process_flipkin_output(kinstr):
    '''Currently just a wrapper for kinemage.process_kinemage().'''
    return kinemage.process_kinemage(kinstr)


def create_object_with_flipkin_coords(mpobj, which_flips='NQ'):
    '''Duplicate the mpobj molecule and apply the flipNQ/H group coordinates.

    PARAMETERS

        mpobj           An MPObject instance

        which_flips     Either 'NQ' or 'H' to specify which flipkin to use

    '''
    flipkin_name = 'flipkin{}'.format(which_flips)
    flipkin = mpobj.kin[flipkin_name]
    flip_group = 'flip{}'.format(which_flips)

    reduced_obj = mpobj.pdb['reduce']
    flipped_obj = mpobj.pdb[flipkin_name]
    cmd.create(flipped_obj, reduced_obj)

    for vl in flipkin.vectorlists():
        # Skip if not coordinates
        if vl[0].vectorlist_name == 'x':
            msg = 'skipping non-coords vectorlist {}'
            logger.debug(msg.format(vl[0].vectorlist_name))
            continue

        # Skip everything except the'flipNQ' (or 'flipH') group
        if vl[0].group[0] != flip_group:
            msg = 'skipping non-{} vectorlist'.format(flip_group)
            logger.debug(msg)
            continue

        logger.debug('begin flipping vectorlist')
        for v in vl:
            if not (v.atom[0] and v.atom[1]):
                msg = 'vector {} was missing atoms: {}'.format(v, v.atom)
                logger.warning(msg)
                continue
            macro = v.macro(1)
            sel = v.sel(1)
            logger.debug('atom to be flipped: {}\n  sel: {}'.format(macro, sel))
            source_coords = v.coords[1]

            target_sel = '{} and {}'.format(flipped_obj, sel)
            msg = 'flipping atom {} to {}'.format(target_sel, source_coords)
            logger.debug(msg)


            ret = cmd.load_coords([source_coords], target_sel)
            if ret == -1:
                success = 0
                a = v.atom[1]

                if a['resn'] == 'HIS':
                    his_h = {
                            'HD1': {
                                'old_h': 'HE2',
                                'old_n': 'NE2',
                                'new_h': 'HD1',
                                'new_n': 'ND1',
                                },
                            'HE2': {
                                'old_h': 'HD1',
                                'old_n': 'ND1',
                                'new_h': 'HE2',
                                'new_n': 'NE2',
                                }
                            }
                    if a['name'] in his_h.keys():
                        # Reduce has switched which N is protonated. Let's move
                        # the old H atom and rename it.
                        resi_sel = '{} and chain {} and resi {}'.format(
                            flipped_obj, a['chain'], a['resi'])
                        h = his_h[a['name']]
                        old_h = '{} and name {}'.format(resi_sel, h['old_h'])
                        old_n = '{} and name {}'.format(resi_sel, h['old_n'])
                        new_h = '{} and name {}'.format(resi_sel, h['new_h'])
                        new_n = '{} and name {}'.format(resi_sel, h['new_n'])

                        # Break the old bond
                        cmd.unbond(old_h, old_n)

                        # Rename the atom
                        cmd.alter(old_h, 'name="{}"'.format(a['name']))

                        # Make the new bond
                        cmd.bond(new_h, new_n)

                        # Retry loading coordinates
                        ret = cmd.load_coords([source_coords], target_sel)
                        if not ret == -1:
                            success = 1

                if not success:
                    msg = 'failed to load coords for {}!'.format(macro)
                    logger.warning(msg)

        logger.debug('end flipping vectorlist')


def flipkin_object(obj):
    '''Run flipkin to generate Asn/Gln and His flip kinemages.

    ARGUMENTS

        obj (str)

            Name of a loaded PyMOL object that has already been passed as an
            argument to `reduce_object()` (`reduce_obj` from the PyMOL command
            line).

    Uses the coordinates output from a previous run of reduce_object().
    '''
    o = get_object(obj)

    # Save a tempfile
    tf = save_to_tempfile(o.reduce_output)

    # Run flipkin to get NQ and H flip kinemages
    flipkinNQ_raw = generate_flipkin_output(tf)
    if not flipkinNQ_raw:
        msg = 'Failed to generate Flipkin NQ output for {}.'.format(obj)
        raise MPException(msg)

    flipkinH_raw = generate_flipkin_output(tf, his=True)
    if not flipkinH_raw:
        msg = 'Failed to generate Flipkin H output for {}.'.format(obj)
        raise MPException(msg)

    # Cleanup
    os.unlink(tf)

    logger.info("Generated Flipkin output for '{}'.".format(obj))

    # Process flipkins
    o.kin['flipkinNQ'] = process_flipkin_output(flipkinNQ_raw)
    o.kin['flipkinH'] = process_flipkin_output(flipkinH_raw)

    create_object_with_flipkin_coords(o, 'NQ')
    create_object_with_flipkin_coords(o, 'H')

    o.draw('flipkinNQ')
    o.draw('flipkinH')

    o.disable_kin('all')
    o.animate()

    logger.info("Stored Flipkin output for '{}'.".format(obj))



###############################################################################
#
#  PROBE
#
###############################################################################

def get_probe_args(pdb_file):
    '''Return arguments for running Probe on the given the PDB filename.'''
    probe_path = 'probe'  # TODO settings.mpgetq('probe_path')
    return [probe_path, '-Quiet', '-Self', 'ALL', pdb_file]


def generate_probe_output(pdbstr):
    '''Generate Probe output from the given PDB string.

    ARGUMENTS
        pdbstr (str)    PDB coordinates from a loaded PyMOL object.

    '''
    # Probe doesn't accept input via STDIN, so we need to write a tempfile.
    tf = save_to_tempfile(pdbstr)
    assert os.path.isfile(tf)
    args = get_probe_args(tf)
    output = run_command(args)
    os.unlink(tf)
    assert not os.path.isfile(tf)
    return output


def process_probe_output(kinstr):
    '''Process Probe output and return lists of dots and clashes.'''

    return kinemage.process_kinemage(kinstr)


def probe_object(obj):
    '''Run Probe on the "Reduce-d" coordinates of a loaded PyMOL object.

    ARGUMENTS

        obj (str)

            Name of a loaded PyMOL object that has already been passed as an
            argument to `reduce_object()` (or `reduce_obj` from the PyMOL
            command line).

    NOTE

        Reduce_object() must be run prior to probe_object() in order to set up
        an MPObject instance in the `objects` dictionary.  Running
        probe_object() on a plain PyMOL object will fail.  Also, accordingly,
        keep in mind that the coordinates probe_object() uses are those of the
        Reduce-modified version.  For an object `myobj`, this will typically
        be `mp_myobj.myobj_reduce`.

    '''


    o = get_object(obj)

    # Clear previous results
    cmd.delete(o.pdb['probe'])  # coords
    cmd.delete(o.get_kin_cgo_group('probe'))  # cgo

    # Create the PDB to use with the 'probe' kinemage
    if o.pdb['userflips'] in cmd.get_names():
        which = 'userflips'
    else:
        which = 'reduce'
    v = cmd.get_view()
    cmd.create(o.pdb['probe'], o.pdb[which])
    cmd.set_view(v)

    pdbstr = o.get_pdbstr('probe')
    output = generate_probe_output(pdbstr)

    # Fail gracefully if probe call returns no output.
    if not output:
        msg = 'Failed to generate Probe output for {}.'.format(obj)
        logger.error(msg)
        return

    logger.info("Generated Probe output for '{}'.".format(obj))

    # Store list of dots and vectors
    o.kin['probe'] = process_probe_output(output)

    # Draw dots and vectors
    o.draw('probe')
    cmd.set_view(v)


# Set up a module-level `objects` storage variable the first time only.
try:  # pragma: no cover
    len(objects)
    logger.info('Using existing MolProbity objects list.')
except:
    objects = {}
    logger.info('Set up MolProbity objects list.')



# ###############################################################################
# #
# #  Set up CLI
# #
# ###############################################################################

cmd.extend('reduce_obj', reduce_object)
cmd.extend('flipkin_obj', flipkin_object)
cmd.extend('probe_obj', probe_object)

logger.info('Finished loading MolProbity plugin.')
