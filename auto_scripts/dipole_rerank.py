#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess
import numpy as np

# Import necessary python modules
import JabberDock as jd
import biobox as bb

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments

# Need number of samples, input style

parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-ir', metavar='receptor_tcl', required=True, help="Initial receptor tcl file output by dock. The name should be consistant between your initial receptor pdb structure and initial receptor dx maps. This is the default behaviour from dock and will only be different if you've moved files")
parser.add_argument('-il', metavar='ligand_tcl', required=True, help="Initial ligand tcl file output by dock. The name should be consistant between your initial ligand pdb structure and initial ligand dx maps. This is the default behaviour from dock and will only be different if you've moved files")
parser.add_argument('-ilpdb', metavar='ligand_pdb', required=True, help="Initial ligand pdb file that corresponds to the initial ligand tcl file input via -il.")
parser.add_argument('-d', metavar='data_file', required=False, default='model_solutions.dat', help="Data file that contains the roto-translation parameters and scores. Default is model_solutions.dat as output by dock.")
parser.add_argument('-o', metavar='new_scores', required=False, default='new_scores.dat', help="Name of re-ranked data file. Note that the output will retain the original order solutions were produced, and will not be ordered according to rank. Default is new_scores.dat")
parser.add_argument('-w', metavar='weight', required=False, default=1.0, help="The weighting applied to the surface complementartiy score when giving the final summed score. The default is 1.0. If you just want the dipole score, set w to 0.0")
parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
args = vars(parser.parse_args())

receptor_tcl = str(args['ir'])
ligand_tcl = str(args['il'])
ligand_pdb = str(args['ilpdb'])
data_file = str(args['d'])
new_scores = str(args['o'])
weight = float(args['w'])

#test to see if files exist
if not os.path.isfile(receptor_tcl):
    raise Exception("ERROR: Receptor tcl file %s not found!"%receptor_tcl)
if not os.path.isfile(ligand_tcl):
    raise Exception("ERROR: Ligand tcl file %s not found!"%ligand_tcl)
if not os.path.isfile(ligand_pdb):
    raise Exception("ERROR: Ligand pdb file %s not found!"%ligand_pdb)
############################################

# create logger
logname = "dipole.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("dipole")
fh = logging.FileHandler(logname)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)

if not args["v"]:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

logger.info("> Loading in files...")

R = jd.Dipole()
L = jd.Dipole()
L_pdb = bb.Molecule()

R.import_dipole(receptor_tcl)
L.import_dipole(ligand_tcl)
L_pdb.import_pdb(ligand_pdb)
data = np.loadtxt(data_file)

COM = L_pdb.get_center()

S = []
no_samples = np.shape(data)[0]

logger.info("> Assessing the dipole alignment...")

for i in range(no_samples):

    logger.info("> Assessing dipoles... %.2f complete"%(float(i) / no_samples))
    
    line = data[i, :]

    # Check to see if dipole alignment is working
    #L_pdb2 = deepcopy(L_pdb)
    #R_mat = jd.geometry.rotate_pdb(L_pdb2, line[3], line[4], line[5], line[6])
    #jd.geometry.translate_pdb(L_pdb2, line[0], line[1], line[2])
    #L_pdb2.write_pdb('test_%i.pdb'%(i))

    L_cp = deepcopy(L)
    L_cp.rotate_dipole(COM, line[3], line[4], line[5], line[6])
    L_cp.translate_dipole(np.asarray((line[0], line[1], line[2])))
    #L_cp.write_dipole('test_%i.tcl'%(i))
    Sd_tmp = R.dipole_score(L_cp)

    if np.isnan(Sd_tmp):
        S.append(weight * line[7])
    else:
        S.append(Sd_tmp + line[7])

logger.info("> Assessment complete. Saving to %s"%(new_scores))
np.savetxt(new_scores, S)

    
