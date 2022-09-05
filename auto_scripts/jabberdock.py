#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess
import numpy as np
import biobox as bb
import JabberDock as jd

# Current path destination:
current_p = str(os.path.dirname(os.path.realpath(__file__)))

# A global script to auto run all the commands to go from a single pdb file and forcefield to a series of predictions. 
# It needs only an input PDB file, as the amber03 forcefield is used by default.
# This script uses all the default settings avaliable for ease of access and use, but keep in mind JabberDock is not foolproof.
# There is some certainty we find the right solution, however it is not concrete (see 2019 paper for more details).

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments
parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-ir', metavar="pdb1", required=True, help='Input PDB file receptor')
parser.add_argument('-il', metavar="pdb2", required=True, help='Input PDB file ligand')
parser.add_argument('-ff', metavar="forcefield", default='amber03', required=False, help='Forcefield reference name (default amber03)')
parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
parser.add_argument('-np', metavar="no_proc", default=1, required=False, help="Number of CPUs to run in parallel. If one, runs in serial (Default 1)")
parser.add_argument('-gpu', metavar="gpu", default=-1, required=False, help="GPU ID to map to for gromacs calculations. 0000 (i.e. first registered one) is the default. But if not specified, then it is assumed you're not using a gpu")
parser.add_argument('-ntomp', metavar="OpemMP_threads", default=0, required=False, help="Number of openMP threads per rank to use (default is zero)")
parser.add_argument('-ff_dat', metavar="Forcefield data file (normally in biobox classes)", default='~/biobox/classes/amber03.dat', required=False, help="Location of dat files used to build the charge maps for STID. There are some examples in biobox/classes, which you may well use as your forcefield. The default is amber03.dat in ~/biobox/classes. This will fail if you've installed biobox elsewhere.")
args = vars(parser.parse_args())

pdb1 = str(args["il"])
pdb2 = str(args["ir"])
ff = str(args["ff"])
n_proc = int(args["np"])
gpu = int(args["gpu"])
ntomp = int(args["ntomp"])
ff_dat = str(args["ff_dat"])
#if variables exist:
if n_proc < 1:
    raise Exception("ERROR: The number of processors you have requested is less than 1!")
elif n_proc > 1 and ntomp < 2:
    raise Exception("ERROR: You have requested multiple MPI threads but no ranks, we recommend an ntomp of 4")

# Test to check if PDB file exists
if not os.path.isfile(pdb1):
    raise Exception("ERROR: PDB file %s not found!"%pdb2)
if not os.path.isfile(pdb2):
    raise Exception("ERROR: PDB file %s not found!"%pdb2)
############################################

# create logger
logname = "jabberdock.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("jabberdock")
fh = logging.FileHandler(logname)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)

if not args["v"]:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

# log initial command
cmd_output_string = sys.argv
cmd_output = ""
for i in range(len(cmd_output_string)):
    cmd_output = cmd_output + cmd_output_string[i] + " "

logger.info(cmd_output)

logger.info("> Beginning gromacs run...")

recep_name = pdb1.split('.')[0]
lig_name = pdb2.split('.')[0]

subprocess.call('python %s/gmx_run.py -i %s -ff %s -gpu %i'%(current_p, pdb1, ff, gpu), shell=True)
subprocess.call('python %s/gmx_run.py -i %s -ff %s -gpu %i'%(current_p, pdb2, ff, gpu), shell=True)

logger.info("> Converting MD output to STID maps...")

subprocess.call('python %s/build_map.py -i %s_sim.pdb -ff %s'%(current_p, recep_name, ff_dat), shell=True)
subprocess.call('python %s/build_map.py -i %s_sim.pdb -ff %s'%(current_p, lig_name, ff_dat), shell=True)

logger.info("> Maps complete, beginning docking procedure...")

subprocess.call('python %s/dock.py -ir %s_sim_map -il %s_sim_map -np %i'%(current_p, recep_name, lig_name, n_proc), shell=True)

#subprocess.call("cd models/", shell=True)

logger.info("> Docking complete, just producing a final file with the scores...")

subprocess.call("python %s/rank.py"%(current_p), shell=True)

top_ten = np.loadtxt('models/ranked_scores.dat')[:10, :]
model_no = top_ten[:, -1]

logger.info("> JabberDock complete. Your top ten models in decending order are %i, %i, %i %i, %i, %i, %i, %i %i and %i"%(model_no[0],model_no[1],model_no[2],model_no[3],model_no[4],model_no[5],model_no[6],model_no[7],model_no[8],model_no[9]))
