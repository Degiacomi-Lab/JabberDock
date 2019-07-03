# JabberDock
JabberDock provides a mechanism to dock two protein STID maps together in conjunction with the POW engine
and BioBox

Note that requires BioBOx as a means to generate the electron density maps, and POW to run the optimisation.

This version of JabberDock ships with a version of biobox intended to allow it to work.

author: Lucas Rudden, l.s.rudden@durham.ac.uk

auto_scripts
Command line features are provided in the auto_scripts folder.
* Help for each one for be provided by typing the command follows by -h, e.g. dock.py -h
* The commands should be accessible in your PYTHONPATH
* More information on each command is provided in the accompanying manual

methods
* In here you will find all the python scripts called by the various commands in auto_scripts
* Each one has documentation written out, this can be read either directly in the file, or via the manual
* There are a few additional commands provided for users which are not used by the commands in auto_scripts...
* For example, the create_pqr_dat in data, which can convert an rtp file into a readable dat pqr file to create STID maps

## INSTALLATION AND REQUIREMENTS ##

JabberDock requires python 2.7 and the following packages (and their associated packages):
* numpy
* biobox
* cython
* scipy
* scikit-image
* subprocess

JabberDock requires the following software to be installed to run:
* gromacs 5.x
* VMD 1.9.x
* POWer

Install with: 'python setup.py install'

## Biobox installation ##
Shipped with this version of JabberDock is a heavily stripped version of biobox containing the functions JabberDock requires.
Please move the biobox folder into whatever directory you like (e.g. your home), just make sure it's in your pythonpath
There is a readme file in the biobox folder containing instructions on requirements and installation instructions.
