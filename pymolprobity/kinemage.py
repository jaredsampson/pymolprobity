'''Kinemage handling for PyMOLProbity plugin.'''

from __future__ import absolute_import
from __future__ import print_function

import copy
import logging
import re

from pymol import cmd

from . import points
from . import utils


logger = logging.getLogger(__name__)



###############################################################################
#
#  KINEMAGE
#
###############################################################################

class Kinemage(object):
    '''Store Kinemage data.'''

    def draw(self, kin_group, dot_mode=0):
        '''Draw the Kinemage dots and clashes.'''
        logger.debug('Drawing Kinemage {}...'.format(self))
        view = cmd.get_view()

        # Set group_auto_mode for easier naming
        gam = cmd.get('group_auto_mode')
        cmd.set('group_auto_mode', 2)

        # Speed up CGO sphere rendering
        if dot_mode == 0:
            cmd.set('cgo_use_shader', 0)
            cmd.set('cgo_sphere_quality', 0)

        # Create dict to store cgolists to be extended
        cgolists = {}

        # Populate a Pointmaster code lookup dict
        pm_lookup = {}
        for pm in self.pointmasters():
            if pm['code'] not in pm_lookup.keys():
                pm_lookup[pm['code']] = pm['label']

        # Get CGO from dotlists
        for dotlist in self.dotlists():
            # Skip non-contact dotlists (if there is such a thing)
            if dotlist[0].dotlist_name != 'x':
                continue

            try:
                dg = dotlist[0].group[0]
            except TypeError:
                dg = 'no_group'  # probe output typically doesn't have a group

            ds = dotlist[0].subgroup[1]
            dm = dotlist[0].master
            for dot in dotlist:
                dpm = pm_lookup[dot.pm]
                dcgo = dot.get_cgo(dot_mode)
                dname = '{}.{}.{}.{}'.format(dg, ds, dm, dpm)
                try:
                    # Extend existing cgo list
                    cgolists[dname].extend(dcgo)
                except KeyError:
                    # Create new cgo list
                    cgolists[dname] = dcgo

        # Get CGO from vectorlists
        # TODO combine this with the dotlist version into a separate function
        for vectorlist in self.vectorlists():
            # Skip non-clash vectorlists (e.g. coordinates)
            if vectorlist[0].vectorlist_name != 'x':
                continue

            try:
                vg = vectorlist[0].group[0]
            except TypeError:
                vg = 'no_group'  # probe output typically doesn't have a group

            vs = vectorlist[0].subgroup[1]
            vm = vectorlist[0].master
            for vector in vectorlist:
                vpm = pm_lookup[vector.pm[0]]  # 2 stored, use first
                vcgo = vector.get_cgo()
                vname = '{}.{}.{}.{}'.format(vg, vs, vm, vpm)
                try:
                    # Extend existing cgo list
                    cgolists[vname].extend(vcgo)
                except KeyError:
                    # Create new cgo list
                    cgolists[vname] = vcgo

        # Create CGO objects
        for name, cgolist in cgolists.items():
            objname = '{}.{}'.format(kin_group, name)
            logger.debug('loading cgo for object {}'.format(objname))
            cmd.load_cgo(cgolist, objname)

        # Restore initial view.
        cmd.set_view(view)

        # Restore initial group_auto_mode setting
        cmd.set('group_auto_mode', gam)

        logger.debug('Finished drawing Kinemage.')

    def get_all_keywords_of_type(self, kw):
        l = []
        for i, k in self.keywords.items():
            if k['keyword'] == kw:
                l.append(k['data'])
        return l

    def get_unique_keywords_of_type(self, kw):
        l = []
        for i, k in self.keywords.items():
            if k['keyword'] == kw:
                if k['data'] not in l:
                    l.append(k['data'])
        return l

    # TODO add filtering to these methods
    # e.g. kin.vectorlists(filter={'group': 'flipNQ'}) to get only those
    # vectorlists in group flipNQ.

    def viewids(self):
        return self.get_all_keywords_of_type('viewid')

    def groups(self):
        return self.get_all_keywords_of_type('group')

    def subgroups(self):
        return self.get_unique_keywords_of_type('subgroup')

    def kin_subgroups(self):
        subgroups = self.subgroups()
        return [sg for sg in subgroups if sg[0] == 'dominant']

    def masters(self):
        return self.get_unique_keywords_of_type('master')

    def pointmasters(self):
        return self.get_unique_keywords_of_type('pointmaster')
        # l = []
        # for i, k in self.keywords.items():
        #     if k['keyword'] == 'pointmaster':
        #         pm = k['data']['label']
        #         if pm not in l:
        #             l.append(pm)
        # return l

    def dotlists(self):
        return self.get_all_keywords_of_type('dotlist')

    def vectorlists(self):
        return self.get_all_keywords_of_type('vectorlist')

    def __init__(self):
        self.keywords = {}



def single_line_keyword_check(lines):
    '''Call this for keywords that should only be a single line.

    If the list has more than one line, print a warning.  If the input is not
    properly formatted, raise a ValueError.
    '''
    if type(lines) is not list:
        msg = 'Expected a list of keyword lines but got: {}'.format(lines)
        raise ValueError(msg)
    if len(lines) > 1:
        msg = '  Only using the first line of multiline keyword: {}'
        kw = lines[0].split()[0]
        logger.warning(msg.format(kw))


VIEWID_RE = re.compile(
        r'''(\d*)viewid '''     # view number
        r'''{([ \*])'''         # is flipped by Reduce? (* or space)
        r'''([\w])'''           # single-letter AA resn
        r'''([ \w]{2,5})'''     # 1-4 digit residue number plus insertion code
        r'''([ \w])'''          # alt
        r'''([ \w])}''')        # chain id
def process_viewid(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = VIEWID_RE.match(line)
    if m:
        return {
                'view_num': (int(m.group(1)) if m.group(1) else 1),
                'flipped': m.group(2) == '*',
                'resn': m.group(3).strip(),
                'resi': m.group(4).strip(),
                'alt': m.group(5).strip(),
                'chain': m.group(6).strip(),
                }
    else:
        logger.warning('Unrecognized viewid format: "{}"'.format(line))
        return None


MASTER_RE = re.compile(r'''master {([^}]*)}''')
def process_master(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = MASTER_RE.match(line)
    if m:
        return utils.slugify(m.group(1))


POINTMASTER_RE = re.compile(
        r'''pointmaster '(\w*)' '''     # pointmaster code
        r'''{([\w\s]*)}'''              # pointmaster label
        r'''(?: (\w+))?''')             # default state
def process_pointmaster(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = POINTMASTER_RE.match(line)
    if m:
        pm = {
                'code': m.group(1),
                'label': utils.slugify(m.group(2)),
                'enable': 0 if m.group(3) == 'off' else 1,  # default to "on"
                }
        return pm


KINEMAGE_RE = re.compile(r'''kinemage ([\w\d]+)''')
def process_kinemage_keyword(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = KINEMAGE_RE.match(line)
    if m:
        return m.group(1)


GROUP_RE = re.compile(r'''group {([^}]*)} (dominant|animate)''')
def process_group(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = GROUP_RE.match(line)
    if m:
        return [m.group(1), m.group(2)]


SUBGROUP_RE = re.compile(r'''subgroup(?: (\w*))? {([^}]*)}(?: (\w+))?(?: (\w+))?''')
def process_subgroup(lines, context):
    single_line_keyword_check(lines)
    line = lines[0]
    m = SUBGROUP_RE.match(line)
    if m:
        return [m.group(1), m.group(2).replace('->', 'to'), m.group(3), m.group(4)]


def process_kinemage(kinstr):
    '''Parse a Probe output string and return a Kinemage object.'''
    KEYWORD_HANDLERS = {
            # MASTERS, ASPECTS, AND COLORS
            'master': process_master,
            'pointmaster': process_pointmaster,

            # KINEMAGES, GROUPS AND SUBGROUPS
            'kinemage': process_kinemage_keyword,
            'group': process_group,
            'subgroup': process_subgroup,

            # LISTS
            'dotlist': points.process_dotlist,
            'vectorlist': points.process_vectorlist,

            # VIEWS
            # may be preceded by a number, e.g. `@2viewid`
            'viewid': process_viewid,
            }

    SKIPPED_KEYWORDS = [
            # METADATA
            'text', 'title', 'copyright', 'caption', 'mage', 'prekin',
            'pdbfile', 'command', 'dimensions', 'dimminmax', 'dimscale',
            'dimoffset',

            # DISPLAY OPTIONS
            'whitebackground', 'onewidth', 'thinline', 'perspective', 'flat',
            'listcolordominant', 'lens',

            # MASTERS, ASPECTS, AND COLORS
            'colorset', 'hsvcolor', 'hsvcolour',

            # LISTS
            'labellist', 'ringlist', 'balllist', 'spherelist', 'trianglelist',
            'ribbonlist', 'marklist', 'arrowlist',

            # VIEWS
            # may be preceded by a number, e.g. `@2span`
            'span', 'zslab', 'center',
            ]

    kin = Kinemage()

    commands = kinstr.lstrip('@').split('\n@')

    # Track context
    context = {
            'kinemage': None,
            'group': None,
            'subgroup': None,
            'animate': 0,
            }
    for i, command in enumerate(commands):
        lines = command.strip().split("\n")
        keyword = lines[0].split(" ")[0]  # First word after "@"
        base_keyword = re.sub(r'\d', '', keyword)  # remove any digits
        if base_keyword in KEYWORD_HANDLERS.keys():

            # Process keyword lines with the function set in KEYWORD_HANDLERS
            logger.debug('Processing keyword {}: {} as {}...'.format(i,
                keyword, base_keyword))
            data = KEYWORD_HANDLERS[base_keyword](lines, copy.copy(context))
            kin.keywords[i] = {
                    'keyword': base_keyword,
                    'data': data
                    }

            logger.debug("  Stored keyword {}: {}.".format(i, keyword))

            # Manage context after kinemage, group, and subgroup keywords
            if base_keyword == 'kinemage':
                context['kinemage'] = data
                msg = 'entering kinemage #{kinemage}'.format(**context)
                logger.debug(msg)

            elif base_keyword == 'group':
                context['group'] = copy.deepcopy(data)
                try:
                    if 'animate' in data:
                        context['animate'] = 1
                    else:
                        context['animate'] = 0
                except TypeError:
                    context['animate'] = 0

                msg = 'entering group: {group} (animate = {animate})'
                logger.debug(msg.format(**context))

            elif base_keyword == 'subgroup':
                context['subgroup'] = copy.deepcopy(data)
                logger.debug('entering subgroup: {subgroup}'.format(**context))


        elif base_keyword in SKIPPED_KEYWORDS:
            logger.debug('Skipping keyword {}: {}'.format(i, keyword))

        else:
            logger.warning('Unknown keyword: {}'.format(keyword))

    return kin







