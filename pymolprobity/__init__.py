'''pymolprobity: A MolProbity plugin for PyMOL

REQUIREMENTS

    The PyMOLProbity Plugin depends on programs from the Richardson lab at Duke
    University, which must be installed locally.  These programs, Reduce,
    Prekin, Flipkin, and Probe, are available for download from the Richardson
    lab website:

        http://kinemage.biochem.duke.edu/software/

    They are available as precompiled binaries for Linux and MacOS (or, in the
    case of Flipkin, a Perl script), and Reduce and Probe are also included
    with Phenix as phenix.reduce and phenix.probe, so if you have an existing
    Phenix installation, you don't need to install them again.

    The easiest way to set the programs up to run is probably to place the
    binaries (or symlinks to them) in a directory on the shell PATH, e.g. in
    `/usr/local/bin`.  Otherwise, save them wherever you like and update the
    appropriate plugin setting from the PyMOL command line:

        mpset reduce_path, /path/to/reduce
        mp_save_settings

    You may also wish to download the Richardson group's "slightly modified
    version of the connectivity table published by the PDB" from the Reduce
    software page above.  This file should be placed in /usr/local.  (If you
    don't, you'll probably get the following error when you run Reduce.

        ERROR CTab(/usr/local/reduce_wwPDB_het_dict.txt): could not open





'''
__author__ = 'Jared Sampson'
__version__ = '0.1'


import logging

import gui


logging.basicConfig(level=logging.INFO)

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'PyMOLProbity',
            label = 'PyMOLProbity',
            command = lambda s=self: gui.PyMOLProbity(s))

