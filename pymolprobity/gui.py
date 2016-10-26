'''GUI for PyMOLProbity plugin.'''
import logging
import Pmw
import Tkinter as tk

from pymol import cmd

import main

logger = logging.getLogger(__name__)


class PyMOLProbity:
    '''The main PyMOLProbity plugin dialog.'''
    def __init__(self, app):

        def get_object_list():
            all_obj = cmd.get_names_of_type('object:molecule')
            # filter to exclude plugin-generated objects
            ret = [o for o in all_obj if not o.startswith('mp_')]
            logger.debug('available objects are: {}'.format(ret))
            return ret

        def update_combo():
            input_combo.setlist(get_object_list())
            input_combo.selectitem(mpobj.get())

        def refresh_page(page):
            try:
                logger.debug('Refreshing page {}.'.format(page))
                # Refresh the combo box in case something has changed
                update_combo(input_combo)
                nb.setnaturalsize()
            except Exception as e:
                logger.debug('Failed to refresh page {}.'.format(page))

        def set_mpobj(val):
            mpobj.set(val)
            cmd.zoom(val, 3)
            logger.debug('mpobj set to {}'.format(val))
            # Get active flipkin
            flipkin = self.flipkin_radio.getvalue()
            o = main.get_object(val)
            if o.kin[flipkin] is not None:
                logger.debug('populating flips group...')
                populate_flips_grp(flipkin)
            else:

                msg = "No {} for object {}...clearing flips.".format(
                        flipkin, val)
                logger.debug(msg)
                clear_child_widgets(flips_grp.interior())

        def reduce_mpobj():
            main.reduce_object(mpobj.get())

        def flipkin_mpobj():
            main.flipkin_object(mpobj.get())
            frame = flipkin_ctl.interior()
            self.flipkin_radio = gen_flipkin_radioselect(frame)
            self.flipkin_group_chk = gen_flipkin_group_checkbuttons(frame)
            self.flipkin_radio.invoke('flipkinNQ')  # always start with NQ
            animate_mpobj()

        def select_flipkin(flipkin):
            o = main.get_object(mpobj.get())
            o.solo_kin(flipkin)

            # Set the proper PDB based on the active group
            try:
                val = self.flipkin_group_chk.getvalue()
            except:
                # first time running, default to NQ
                val = 'flipkinNQ'
            logger.debug('self.flipkin_group_chk = {} during select_flipkin'.format(val))
            if 'flip' in val:
                o.solo_pdb(flipkin)
                if 'reduce' in val:
                    o.enable_pdb('reduce')
            else:
                o.solo_pdb('reduce')

            # Set up the view controls
            frame = flipkin_ctl.interior()
            self.flipkin_subgroups_chk = gen_flipkin_subgroup_checkbuttons(frame)
            self.flipkin_masters_chk = gen_flipkin_master_checkbuttons(frame)
            # self.flipkin_pointmasters_chk = gen_flipkin_pointmaster_checkbuttons(frame)

            # Show and populate the flips group for this flipkin
            populate_flips_grp(flipkin)
            nb.setnaturalsize()

        def enable_flipkin_group(group, enable=True):
            o = main.get_object(mpobj.get())

            # use full names (e.g. 'flipkinNQ') for PDBs
            if group.startswith('flip'):
                pdb = self.flipkin_radio.getvalue()
                suffix = pdb.replace('flipkin','')  # yields 'NQ' or 'H'
                assert suffix in ['NQ', 'H']
                kin_grp = '{}{}'.format(group, suffix)
            else:
                # unless something is wrong with the flipkin, group is 'reduce'
                pdb = group
                kin_grp = group

            # enable (or disable if enable=False)
            if enable:
                logger.debug('enabling {} and {}'.format(kin_grp, pdb))
                o.enable_flipkin_group(kin_grp)
                o.enable_pdb(pdb)
            else:
                logger.debug('disabling {} and {}'.format(kin_grp, pdb))
                o.disable_flipkin_group(kin_grp)
                o.disable_pdb(pdb)

        def animate_mpobj():
            o = main.get_object(mpobj.get())
            # Ensure only one is checked
            val = self.flipkin_group_chk.getvalue()
            if len(val) == 1:
                # Toggle both buttons
                self.flipkin_group_chk.invoke('reduce')
                self.flipkin_group_chk.invoke('flip')
            elif len(val) in [0, 2]:
                # If they're both checked or unchecked, only toggle one
                self.flipkin_group_chk.invoke('reduce')

        def toggle_subgroup(subgroup, enable):
            SEL = {
                    '1->2': 'mp_*.1to2.*',
                    '2->1': 'mp_*.2to1.*',
                    }
            if enable:
                cmd.enable(SEL[subgroup])
            else:
                cmd.disable(SEL[subgroup])

        def toggle_master(master, enable):
            SEL = {
                    'vdw contact': 'mp_*vdw_contact',
                    'small overlap': 'mp_*small_overlap',
                    'bad overlap': 'mp_*bad_overlap',
                    'H-bonds': 'mp_*H_bonds',
                    }
            if enable:
                cmd.enable(SEL[master])
                # Keep buttons on Flipkin and Probe tabs in sync.
                try:
                    self.flipkin_masters_chk.button(master).select()
                    self.probe_masters_chk.button(master).select()
                except:
                    # don't choke if the button doesn't exist
                    pass

            else:
                cmd.disable(SEL[master])
                try:
                    self.flipkin_masters_chk.button(master).deselect()
                    self.probe_masters_chk.button(master).deselect()
                except:
                    pass

        def toggle_pointmaster(pointmaster, enable):
            pass

        def zoom_and_refresh(sel):
            cmd.orient(sel)
            cmd.zoom(sel, 3)
            cmd.refresh()

        def clear_child_widgets(parent):
            # TODO may be better to keep the widgets and use grid_forget()
            for widget in parent.winfo_children():
                widget.destroy()
            logger.debug('cleared contents of {}'.format(parent))

        def populate_flips_grp(flipkin):
            parent = flips_grp.interior()
            # Clear any existing child widgets
            clear_child_widgets(parent)
            # Recreate the scrolling frame to house the flip rows
            frame = gen_scrolling_frame(parent)

            # Column headers
            gen_flips_list_header(frame)

            # Repopulate with flips
            o = main.get_object(mpobj.get())
            if not o.views[flipkin]:
                kin = o.kin[flipkin]
                o.views[flipkin] = kin.viewids()
            logger.debug('begin populating flip rows')
            for i, v in enumerate(o.views[flipkin]):
                gen_flips_list_row(frame, o, flipkin, i+1, v)
            logger.debug('finished populating flip rows')
            nb.setnaturalsize()

        def save_flip_selections():
            '''Save the user-selected flips to a new PDB object.'''
            o = main.get_object(mpobj.get())
            flipkin = self.flipkin_radio.getvalue()
            reduce_obj = o.pdb['reduce']
            flipkin_obj = o.pdb[flipkin]
            userflips_obj = o.pdb['userflips']
            # Create a new userflips object even if it already exists in case
            # previous selections have been changed.
            v = cmd.get_view()
            cmd.create(userflips_obj, reduce_obj)
            cmd.set_view(v)

            for i, v in enumerate(o.views[flipkin]):
                # If reduce value and user value are different, copy the
                # coordinates from the current flipkin molecule
                if v['reduce_chk_val'].get() != v['user_chk_val'].get():
                    # Do it the hard way, by combining objects.  This is plenty
                    # fast (we typically won't have too many flips to switch)
                    # and doesn't result in atom name mismatch errors for
                    # differently protonated HIS residues the way the
                    # load_coords method does.
                    flipped_sel = '({} and chain {} and resi {})'.format(
                            flipkin_obj, v['chain'], v['resi'])
                    userflips_sel = '({} and not (chain {} and resi {}))'.format(
                            userflips_obj, v['chain'], v['resi'])
                    combined_sel = '{} or {}'.format(
                            userflips_sel, flipped_sel)
                    v = cmd.get_view()
                    cmd.create(userflips_obj, combined_sel)
                    cmd.set_view(v)
                    msg = 'added flip for {} to {}'.format(flipped_sel,
                            userflips_obj)
                    logger.debug(msg)
            o.solo_pdb('userflips')

        def probe_mpobj():
            main.probe_object(mpobj.get())
            frame = probe_ctl.interior()
            # Set up the view controls
            # self.probe_group_chk = gen_flipkin_group_checkbuttons(frame)
            # self.probe_subgroups_chk = gen_probe_subgroup_checkbuttons(frame)
            self.probe_masters_chk = gen_probe_master_checkbuttons(frame)
            # self.probe_pointmasters_chk = gen_probe_pointmaster_checkbuttons(frame)
            nb.setnaturalsize()



        ###############################################################################
        #
        #  Widget Generators
        #
        ###############################################################################

        def gen_input_obj_combo(parent):
            '''Show the available loaded PyMOL objects in a ComboBox.'''
            obj_list = get_object_list()

            # Stop if there are no available objects.
            if len(obj_list) == 0:
                w = tk.Label(parent, text='No objects loaded!')
                w.grid(sticky='ew')
                return None

            # Create a frame container
            frame = tk.Frame(parent)
            frame.pack(fill='x')

            # Create and populate the combo box
            combo = Pmw.ComboBox(frame,
                    label_text = 'Base object name:',
                    labelpos = 'w',
                    selectioncommand = set_mpobj,
                    scrolledlist_items = obj_list,
                    )
            combo.grid(row=0, sticky='ew')

            # Create a refresh button
            btn = tk.Button(frame, text='Refresh List', command=update_combo)
            btn.grid(row=0, column=1)

            # Select first object by default if not previously set
            if not mpobj.get():
                combo.selectitem(obj_list[0])
                mpobj.set(obj_list[0])
            return combo

        def gen_flipkin_radioselect(parent):
            # destroy the previous one if it exists
            try:
                self.flipkin_radio.destroy()
            except AttributeError:
                pass

            # Radio buttons to select a flipkin.
            select = Pmw.RadioSelect(parent,
                    labelpos='w',
                    label_text='Select a flipkin:',
                    buttontype='radiobutton',
                    command=select_flipkin)
            select.add('flipkinNQ')
            select.add('flipkinH')
            select.grid(sticky='ew')
            return select

        def gen_flipkin_group_checkbuttons(parent):
            # destroy the previous version, if it exists
            try:
                self.flipkin_group_chk_frame.destroy()
            except AttributeError:
                pass
            frame = tk.Frame(parent)
            frame.grid(sticky='ew')
            self.flipkin_group_chk_frame = frame

            # Checkbuttons to choose flipkin group(s).
            select = Pmw.RadioSelect(frame,
                    labelpos='w',
                    label_text='Select flipkin group(s):',
                    buttontype='checkbutton',
                    command=enable_flipkin_group)
            select.add('reduce')
            select.add('flip')
            select.grid(row=0, column=1)
            # Animate button
            btn = tk.Button(frame, text='Animate', command=animate_mpobj)
            btn.grid(row=0, column=2)
            # Return only the RadioSelect widget for access to values & names
            return select

        def gen_flipkin_subgroup_checkbuttons(parent):
            # destroy the previous version
            try:
                self.flipkin_subgroups_chk.destroy()
            except AttributeError:
                pass

            chk = Pmw.RadioSelect(parent,
                    buttontype='checkbutton',
                    labelpos='w',
                    label_text='Subgroups:',
                    command=toggle_subgroup)
            chk.grid(sticky='ew')
            for btn in ['1->2', '2->1']:
                chk.add(btn)
                chk.invoke(btn)
            # o = main.get_object(mpobj.get())
            # kinNQ = o.kin['flipkinNQ']
            # kinNQ_sg = [x[1] for x in kinNQ.kin_subgroups()]
            # kinH = o.kin['flipkinH']
            # kinH_sg = [x[1] for x in kinH.kin_subgroups()]

            # Get sorted list union of all subgroups (should both be the same)
            # subgroups = sorted(set().union(kinNQ_sg, kinH_sg))
            # for sg in subgroups:
            #     chk.add(sg)
            return chk

        def gen_flipkin_master_checkbuttons(parent):
            # destroy the previous version if it exists
            try:
                self.flipkin_masters_chk.destroy()
            except AttributeError:
                pass

            chk = Pmw.RadioSelect(parent,
                    buttontype='checkbutton',
                    labelpos='w',
                    label_text='Masters:',
                    command=toggle_master)
            chk.grid(sticky='w')

            # Bonds
            # bond_masters = ['mainchain', 'sidechain', 'hydrogens', 'hets',
            #         'waters']

            # Contacts
            contact_masters = ['vdw contact', 'small overlap',
                    'bad overlap', 'H-bonds']
            for btn in contact_masters:
                chk.add(btn)
                chk.invoke(btn)

            return chk

        # def gen_flipkin_pointmaster_checkbuttons(parent):
        #     # destroy the previous version if it exists
        #     try:
        #         self.flipkin_pointmasters_chk.destroy()
        #     except AttributeError:
        #         pass

        #     chk = Pmw.RadioSelect(parent,
        #             buttontype='checkbutton',
        #             labelpos='w',
        #             label_text='Pointmasters:',
        #             command=toggle_pointmaster)
        #     chk.grid(sticky='ew')
        #     o = main.get_object(mpobj.get())
        #     kinNQ = o.kin['flipkinNQ']
        #     kinNQ_sg = [x[1] for x in kinNQ.pointmasters()]
        #     kinH = o.kin['flipkinH']
        #     kinH_sg = [x[1] for x in kinH.pointmasters()]

        #     # Get sorted list union of all subgroups (should both be the same)
        #     masters = list(set().union(kinNQ_sg, kinH_sg))
        #     for m in masters:
        #         chk.add(m)
        #     return chk

        # Set up scrollbar, canvas and data frame
        # Based on http://stackoverflow.com/a/3092341
        def gen_scrolling_frame(parent, **canvas_options):
            '''Return a frame + scrollbar packed within a Tkinter Canvas.'''
            def onFrameConfigure(canvas):
                '''Reset the scroll region to encompass the inner frame.'''
                canvas.configure(scrollregion=canvas.bbox('all'))
            canvas = tk.Canvas(parent, **canvas_options)
            frame = tk.Frame(canvas)
            vsb = tk.Scrollbar(parent, orient='vertical', command=canvas.yview)
            canvas.configure(yscrollcommand=vsb.set)
            vsb.pack(side='right', fill='y')
            canvas.pack(side='left', fill='both', expand=1)
            canvas.create_window((0,0), window=frame, anchor='nw')
            frame.bind('<Configure>',
                    lambda event, canvas=canvas: onFrameConfigure(canvas))
            return frame

        def gen_flips_list_header(frame):
            # Formatting for all labels
            labels = ['Flipped Residue', 'Zoom', 'Flipped by Reduce',
                    'User Selection']  # TODO add 'Current'?
            for i, label in enumerate(labels):
                w = tk.Label(frame, text=label,
                        relief='ridge', fg='#ffffff', bg='#555555')
                w.grid(row=0, column=i,
                        sticky='nsew', ipadx=4, ipady=4)

        def gen_flips_list_row(frame, o, flipkin, i, v):
            '''Generate a row in the flips list frame for the given view.

            PARAMETERS

                frame   The target frame where the widgets should be drawn.

                o       The MPObject instance to which the viewid belongs.

                i       The row index where the widgets should be gridded.

                v       The viewid dictionary.

            '''
            # TODO make this conversion when processing flipkin
            RESN = {'N': 'ASN', 'Q': 'GLN', 'H': 'HIS'}
            resn = RESN[v['resn']]
            macro = '{}/{}`{}'.format(v['chain'], resn, v['resi'])
            if v['flipped']:
                label_text = '{}*'.format(macro)
            else:
                label_text = macro

            # Flip macro label
            label = tk.Label(frame, text=label_text)
            label.grid(row=i, column=0)
            logger.debug('made label', label_text)

            # Zoom button
            zoom_obj = o.pdb[flipkin]
            logger.debug('got zoom_obj', zoom_obj)
            zoom_sel = '/{}//{} and sidechain'.format(zoom_obj, macro)
            logger.debug('got zoom_sel', zoom_sel)
            zoom_cmd = (lambda s=zoom_sel: zoom_and_refresh(s))
            logger.debug('got zoom_cmd', zoom_cmd)
            zoom_btn = tk.Button(frame, text='Zoom', command=zoom_cmd)
            zoom_btn.grid(row=i, column=1)

            # Flipped by Reduce checkbutton
            # Store value in a Tkinter IntVar (first time only)
            if 'reduce_chk_val' not in v.keys():
                v['reduce_chk_val'] = tk.IntVar()
                v['reduce_chk_val'].set(v['flipped'])
            reduce_chk = tk.Checkbutton(frame, variable=v['reduce_chk_val'],
                    state='disabled')
            reduce_chk.grid(row=i, column=2)

            # User-selected flip checkbutton
            # Create user flip selection storage variable, default to same
            # as 'reduce' value.
            if 'user_chk_val' not in v.keys():
                v['user_chk_val'] = tk.IntVar()
                v['user_chk_val'].set(v['flipped'])
            user_chk = tk.Checkbutton(frame, variable=v['user_chk_val'])
            user_chk.grid(row=i, column=3)

            logger.debug('...added flip {}'.format(macro))

        def gen_probe_master_checkbuttons(parent):
            # destroy the previous version if it exists
            try:
                self.probe_masters_chk.destroy()
            except AttributeError:
                pass

            chk = Pmw.RadioSelect(parent,
                    buttontype='checkbutton',
                    labelpos='w',
                    label_text='Masters:',
                    command=toggle_master)
            chk.grid(sticky='w')

            # Contacts
            contact_masters = ['vdw contact', 'small overlap',
                    'bad overlap', 'H-bonds']
            for btn in contact_masters:
                chk.add(btn)
                chk.invoke(btn)

            return chk


        #######################################################################
        #
        #  Instance Variable
        #
        #######################################################################

        # The current MPObject name
        mpobj = tk.StringVar()


        #######################################################################
        #
        #  Main dialog & tabs
        #
        #######################################################################

        # Main window
        root = app.root
        dialog = Pmw.Dialog(root, title="PyMOLProbity", buttons=[])

        # Add the mpobj selection combo box
        input_combo = gen_input_obj_combo(dialog.interior())

        # Add notebook with tabs
        nb = Pmw.NoteBook(dialog.interior(),
                raisecommand=refresh_page)
        nb.pack(fill='both', expand=1)

        # Page names
        reduce_page_name = 'Add Hydrogens'
        flipkin_page_name = 'Review Flips'
        probe_page_name = 'Visualize Contacts'

        # Create pages
        reduce_page = nb.add(reduce_page_name)
        flipkin_page = nb.add(flipkin_page_name)
        probe_page = nb.add(probe_page_name)


        #######################################################################
        #
        #  Reduce Page
        #
        #######################################################################

        reduce_button = tk.Button(reduce_page, text='Run Reduce',
                command=reduce_mpobj)
        reduce_button.pack(anchor='n', fill='x', expand=1)


        #######################################################################
        #
        #  Flipkin Page
        #
        #######################################################################

        flipkin_button = tk.Button(flipkin_page, text='Run Flipkin',
                command=flipkin_mpobj)
        flipkin_button.pack(anchor='n', fill='x', expand=1)

        # Controls group
        flipkin_ctl = Pmw.Group(flipkin_page, tag_text='View Options')
        flipkin_ctl.pack(anchor='n', fill='x', expand=1)
        self.flipkin_radio = None
        self.flipkin_group_chk = None

        # Flips group
        flips_grp = Pmw.Group(flipkin_page, tag_text='Select Flips')
        flips_grp.pack(anchor='n', fill='x', expand=1)

        # Save flip selections
        save_flips_btn = tk.Button(flipkin_page,
                text='Save Selections', command=save_flip_selections)
        save_flips_btn.pack(anchor='n', fill='x', expand=1)


        #######################################################################
        #
        #  Probe Page
        #
        #######################################################################

        probe_button = tk.Button(probe_page, text='Run Probe',
                command=probe_mpobj)
        probe_button.pack(anchor='n', fill='x', expand=1)

        # Controls group
        probe_ctl = Pmw.Group(probe_page, tag_text='View Options')
        probe_ctl.pack(anchor='n', fill='x', expand=1)
        self.probe_radio = None
        self.probe_group_chk = None
