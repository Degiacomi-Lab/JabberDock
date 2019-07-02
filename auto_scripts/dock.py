#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess
import datetime
now = datetime.datetime.now()
import numpy as np

# Current path destination:
current_p = str(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(current_p)
import build_scripts as bs

# Import necessary python modules
import JabberDock as jd

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments
parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-ir', metavar='receptor name', required=True, help="Name of the receptor. The name should be consistant between your singular pdb structure based around your STID maps. This is the default behaviour from build_map and will only be different if you've moved files")
parser.add_argument('-il', metavar='ligand name', required=True, help="Name of the ligand. The name should be consistant, as above")
parser.add_argument('-iso', metavar='isovalue', required=False, default=0.43, help="Isovalue to use for surface complementartiy of STID map (default is the benchmarked 0.43)")
parser.add_argument('-d', metavar='dist', required=False, default=1.6, help="Distances between corresponding surfaces to consider when using the surface complementarity scoring function")
parser.add_argument('-l', metavar='logfile', required=False, default='pow_log.dat', help="Name of log file to keep POW information")
parser.add_argument('-ns', metavar='no_samples', required=False, default=300, help="Number of docked predicted results you want to be returned")
parser.add_argument('-b', metavar='buffer', required=False, default = 0., help="Buffer region to use when exploring the cartisian conformational space around the receptor (default is 0. ang.)")
parser.add_argument('-m', metavar='space_multiplier', required=False, default= 1.5, help="Multiplier to searching the conformational space. So the default of 2 means twice the size of the receptor is allowed to be explored by the ligand.")
parser.add_argument('-a', metavar='angle', required=False, default = 180., help="Angle to twist molecule around (default / recommended unless you want to restrict the exploration space, 180 deg.)")
parser.add_argument('-n', metavar='name', required=False, default='model_solutions', help="filename for the .dat file with the roto-translations and scores at the end (default is model_solutions)")
parser.add_argument('-r', metavar='restart', required=False, default=0, help="If this is a restart run (i.e. you're continuing from a point in the docking simulation). The default is 0 (i.e. False). If you are continuing, set -r 1")
parser.add_argument('-np', metavar='no_processors', required=False, default=1, help="How many processes to use if running in MPI, default is 1, i.e. serial")
parser.add_argument('-tcl', metavar='dip_map', required=False, default=1, help="Do you want the initial docking position of the tcl dipole map along with its geometric equivilents in the dx and pdb maps. If you want to do any dipole post processing, the answer should be 1 (as is default). The tcl file should have the same name as your STID maps pdb file, e.g. receptor.tcl")

parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
args = vars(parser.parse_args())
receptor = str(args["ir"])
ligand = str(args["il"])
iso = float(args["iso"])
dist = float(args["d"])
logfile = str(args["l"])
no_samples = int(args["ns"])
buff = float(args["b"])
times = float(args["m"])
angle = float(args["a"])
no_proc = int(args["np"])
name = str(args["n"])
if int(args["tcl"]) == 1:
    tcl = True
else:
    tcl = False
if int(args["r"]) == 1:
    restart = True
else:
    restart = False

##################### !!!!!! #################
# Change setting here if POW.py is located in a different place than your user home folder
# It also assumes that the POW installed is the most recent (v2)
pow_loc = "~/POW_v2"
####################

if no_proc == 1:
    mpi = False
else:
    mpi = True

# create logger
logname = "dock.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("dock")
fh = logging.FileHandler(logname)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)

if not args["v"]:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

logger.info("> Setting up input script and defining boundary conditions...")
boundary = jd.geometry.get_minmax_crd(receptor + '.pdb', fname_times = times, fname_buff = buff)

if restart:
    _ = bs.powrun(boundary[0], boundary[1], boundary[2], receptor, ligand, iso=iso, dist = dist, log=logfile, angle=angle, no_samples=no_samples, file_name = name, tcl =tcl, restart = True)
else:
    _ = bs.powrun(boundary[0], boundary[1], boundary[2], receptor, ligand, iso=iso, dist = dist, log=logfile, angle=angle, no_samples=no_samples, file_name = name, tcl =tcl)    

logger.info("> Beginning POW run")

# Creating a models folder to put found solutions if none exists and this isn't a restart run
if not restart:
    if not os.path.isdir('models'):
        subprocess.call('mkdir models', shell=True)
    else:
        logger.info('> Found a models folder, renaming to models_old_%d_%d_%d'%(now.day,now.month,now.year))
        subprocess.call('mv models models_old_%d_%d_%d'%(now.day,now.month,now.year), shell=True)
        subprocess.call('mkdir models', shell=True)

if mpi:
    subprocess.call("mpirun -np %i %s/POW.py input_ensemble > time_log.dat"%(no_proc, pow_loc), shell=True)
else:
    subprocess.call("%s/POW.py input_ensemble > time_log.dat"%(pow_loc), shell=True)

logger.info("%i ligand models produced by POW in folder models. The data file with details on roto-translations and scores is held in %s.dat.\nThe results are returned unranked, but each model number relates to the corresponding line in the data file.\nAn additional dipole re-ranking procedure can now be invoked, or you can rank the solutions as is.\nNote that the solutions are geometrically oriented relative to the initial_receptor file (i.e. the starting points were moved to the origin)"%(no_samples, name))

