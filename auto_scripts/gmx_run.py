#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess

import numpy as np

# Current path destination:
current_p = str(os.path.dirname(os.path.realpath(__file__)))

# Import necessary python modules
sys.path.append(current_p)
import build_scripts as bs

""" Gromacs run scipt"""

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments
parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-i', metavar="pdb", required=True, help='Input PDB file')
parser.add_argument('-ff', metavar="ff", default='amber03', required=False, help='Forcefield reference name (default amber03)')
parser.add_argument('-ts', metavar="timestep", default=2, required=False, help='Timestep of simulation (default 2fs)')
parser.add_argument('-ns', metavar="no_steps", default=300000, required=False, help='Number of steps in MD prod. sim. (default 300000 == 600 ps)')
parser.add_argument('-dt', metavar="dumptime", default=2500, required=False, help='How many steps to dump a frame (default 5 ps)')
parser.add_argument('-t', metavar="temp", default=310.15, required=False, help='Temperate of simulation (default 310.15K)')
parser.add_argument('-p', metavar="press", default=1.0, required=False, help='Pressure of simulation (default 1 bar)')
parser.add_argument('-s', metavar="skip", default=1, required=False, help='How many frames to skip during parsing (default None)')
parser.add_argument('-np', metavar="no_proc", default=1, required=False, help="Number of CPUs to run in parallel. If set to 1 or not specified, use GROMACS default behavior (default is 1)")
parser.add_argument('-gpu', metavar="gpu", default=-1, required=False, help="GPU ID to map to. 0000 (i.e. first registered one) is the default. But if not specified, then it is assumed you're not using a gpu")
parser.add_argument('-minim', metavar="min_script", default=current_p + "/minim.mdp", required=False, help="Path to minimisation script, default is the minim script in the auto_scripts folder")
parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
parser.add_argument('-ntomp', metavar="OpemMP_threads", default=0, required=False, help="Number of openMP threads per rank to use (default is zero)")
args = vars(parser.parse_args())

pdb = str(args["i"])
ff = str(args["ff"])
timestep = int(args["ts"])
no_steps = int(args["ns"])
dump_time = int(args["dt"])
temp = float(args["t"])
press = float(args["p"])
skip = int(args["s"])
n_proc = int(args["np"])
gpu = int(args["gpu"])
min_script = str(args["minim"])
ntomp = int(args["ntomp"])
ion_script = str(current_p + '/ions.mdp')

# If variables exist:
if n_proc != 1:
    gmx_cmd = '/usr/local/gromacs/bin/gmx_mpi'
elif n_proc < 1:
    raise Exception("ERROR: The number of processors you have requested is less than 1!")
else:
    gmx_cmd = '/usr/local/gromacs/bin/gmx'

if n_proc > 1 and ntomp < 2:
    raise Exception("ERROR: You have requested multiple MPI threads but not enough openMP threads.")

if args["gpu"] != -1:
    gpu = str(args["gpu"])
    gpu_cmd = '-gpu_id ' + gpu
else:
    gpu_cmd = ''

# Test to check if PDB file exists
if not os.path.isfile(pdb):
    raise Exception("ERROR: PDB file %s not found!"%pdb)

############################################

# create logger
logname = "gmx_run.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("gmx_run")
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

logger.info("> Initiating MD setup...")

# Verbose or not?
if not args["v"]:
    quiet = '-quiet'
else:
    quiet = ''

# Create gro files
subprocess.call("%s %s pdb2gmx -f %s -o filename.gro -ff %s -water tip3p -ignh"%(gmx_cmd, quiet, pdb, ff), shell=True)

# Get neutralisation value for system
neu_charge = int(subprocess.check_output(current_p + "/neutralise.sh", shell=True))

# Solvate and neutralise
subprocess.call("%s %s editconf -f filename.gro -o filename_proc.gro -c -d 1.0 -bt cubic"%(gmx_cmd, quiet), shell=True)
subprocess.call("%s %s solvate -cp filename_proc.gro -cs spc216.gro -o filename_solv.gro -p topol.top"%(gmx_cmd, quiet), shell=True)
subprocess.call("%s %s grompp -f %s -c filename_solv.gro -p topol.top -o ions.tpr -maxwarn 1"%(gmx_cmd, quiet, ion_script), shell=True)
    
if neu_charge > 0:
    subprocess.call("echo 'SOL' | %s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -nname CL -nn %i"%(gmx_cmd, quiet, neu_charge), shell=True)
elif neu_charge < 0:
    subprocess.call("echo 'SOL' | %s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -np %i -nname CL"%(gmx_cmd, quiet, -1 * neu_charge), shell=True)
else:
    subprocess.call("%s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -nname CL"%(gmx_cmd, quiet), shell=True)

# Begin EM
logger.info("> Beginning energy minimisation")

subprocess.call("%s %s grompp -f %s -c filename_ions.gro -p topol.top -o em.tpr"%(gmx_cmd, quiet, min_script), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -deffnm em"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -deffnm em"%(gmx_cmd, quiet), shell=True)

# Begin Equib
logger.info("> Minimisation complete. Initialising equilibriation...")
_ = bs.nvtrun(temp)

subprocess.call("%s %s grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr"%(gmx_cmd, quiet), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -s nvt.tpr"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -s nvt.tpr"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv confout.gro nvt.gro", shell=True)

# Begin production
logger.info("> Equilibiration complete. Initialising production...")
_ = bs.nptrun(timestep, no_steps, dump_time, temp, press)

subprocess.call("%s %s grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr"%(gmx_cmd, quiet), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -s npt.tpr"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -s npt.tpr"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv traj.trr npt.trr", shell=True)

# Parsing phase
logger.info("> Gromacs run successfully completed, beginning parsing the data")
logger.info("> Converting trajectory to multipdb")

subprocess.call("echo '1' | %s %s trjconv -f npt.trr -s npt.tpr -skip %i -o parse_traj.trr -pbc mol -ur compact"%(gmx_cmd, quiet, skip), shell=True)

# Run the tcl script to produce the multipdb
vmd_file = current_p + "/no_wat.tcl"
subprocess.call("vmd -dispdev text -e %s"%(vmd_file), shell=True)
name_tmp = pdb.split('.')[0]
subprocess.call("mv sim.pdb %s_sim.pdb"%(name_tmp), shell=True)
logger.info("> Gromacs run completed successfully. Check %s_sim.pdb for the final result and input to the next stage"%(name_tmp))
