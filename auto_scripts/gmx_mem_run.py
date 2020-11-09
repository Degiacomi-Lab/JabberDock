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
parser.add_argument('-i', metavar="pdb", default='fixed_trim.pdb', required=False, help='Input PDB file (default fixed_trim.pdb following create_mem_system.py)')
parser.add_argument('-ff', metavar="ff", default='amber03_slipid', required=False, help='Forcefield reference name (default amber03_slipid)')
parser.add_argument('-ts', metavar="timestep", default=2, required=False, help='Timestep of simulation in fs (default 2fs)')
parser.add_argument('-lip', metavar="lipid", default="POPE", required=False, help='Membrane compostiion is (default POPE). NOTE, this only affects the GROMACS groupings, we recommend using POPE regardless, but POPE is the naming convention of slipid. If you use charmm etc., you will need to make sure the lipid is name is correct in -ff. ')
parser.add_argument('-ns', metavar="no_steps", default=5000000, required=False, help='Number of steps in MD prod. sim. (default 5000000 == 10000 ps)')
parser.add_argument('-dt', metavar="dumptime", default=2500, required=False, help='How many steps to dump a frame (default 5 ps)')
parser.add_argument('-t', metavar="temp", default=310.15, required=False, help='Temperature of simulation (default 310.15K)')
parser.add_argument('-p', metavar="press", default=1.0, required=False, help='Pressure of simulation (default 1 bar)')
parser.add_argument('-s', metavar="skip", default=2, required=False, help='How many frames to skip during parsing (default every 2 (10 ps))')
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
lipid = str(args["lip"])
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
logname = "gmx_mem_run.log"
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

# setup condition that an amber script with slipid is being used (for hacking calls)
if "amber" in ff.lower() and "slipid" in ff.lower():
    amber_con = True
else:
    amber_con = False

logger.info("> Initiating MD setup...")

# Verbose or not?
if not args["v"]:
    quiet = '-quiet'
else:
    quiet = ''

# Create gro files
logger.info("> Processing PDB files")
subprocess.call("%s %s pdb2gmx -f %s -o filename.gro -ff %s -water tip3p"%(gmx_cmd, quiet, pdb, ff), shell=True)

# hack topol if amber and slipid
if amber_con:
    logger.info("> Hacking topology to accomodate slipid usage with amber...")
    subprocess.call(current_p + "/hack_topol.sh", shell=True)

# Get neutralisation value for system
neu_charge = int(subprocess.check_output(current_p + "/neutralise.sh", shell=True))

# Neutralise pre solvated system
logger.info("> Creating neutralised system...")
subprocess.call("%s %s editconf -f filename.gro -o filename_proc.gro -c -d 0.0 -bt cubic"%(gmx_cmd, quiet), shell=True)

subprocess.call("%s %s grompp -f %s -c filename_proc.gro -p topol.top -o ions.tpr -maxwarn 1"%(gmx_cmd, quiet, ion_script), shell=True)

# record if the ion grouping scripts are necessary (i.e. does the system need balencing?) - different mdp files if do
if neu_charge > 0:
    ions_present = True
    subprocess.call("echo 'SOL' | %s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -nname CL -nn %i"%(gmx_cmd, quiet, neu_charge), shell=True)
elif neu_charge < 0:
    ions_present = True
    subprocess.call("echo 'SOL' | %s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -np %i -nname CL"%(gmx_cmd, quiet, -1 * neu_charge), shell=True)
else:
    ions_present = False
    subprocess.call("%s %s genion -s ions.tpr -o filename_ions.gro -p topol.top -pname NA -nname CL"%(gmx_cmd, quiet), shell=True)

# Begin EM
logger.info("> Beginning energy minimisation")

subprocess.call("%s %s grompp -f %s -c filename_ions.gro -p topol.top -o em.tpr"%(gmx_cmd, quiet, min_script), shell=True)

if n_proc != 1:
    ############### PLEASE NOTE BELOW ############
    # an alternative to this is "gmx mdrun ntmpi %i -ntomp %i %s -s em.tpr" (or equivilent) (etc., depending on your GROMACS installation)
    # e.g. gromacs with GPU installation will usually need the above command (amend as necessary to all 4 calls to mdrun).
    ##############################################
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -deffnm em"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -deffnm em"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv confout.gro em.gro", shell=True)

# Begin Equib
logger.info("> Minimisation complete. Initialising equilibriation...")
logger.info("> Occasionally this fails with LINCS errors from lipid clashes. If this occurs, rerun create_system.py (delete bilayer_*.pdb file first)")
if ions_present:
    _ = bs.nptrun_mem(timestep=2, no_steps=10000000, dump_time=10000, ref_temp=temp, ref_press=press, continuation="no", constrain=True, lipid=lipid)
else:
    _ = bs.nptrun_mem(timestep=2, no_steps=10000000, dump_time=10000, ref_temp=temp, ref_press=press, continuation="no", constrain=True, ions=False, lipid=lipid)

subprocess.call("%s %s grompp -f npt_constrain.mdp -c em.gro -r em.gro -p topol.top -o npt_constrain.tpr"%(gmx_cmd, quiet), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -s npt_constrain.tpr"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -s npt_constrain.tpr"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv confout.gro npt_constrain.gro", shell=True)

# second stage (no constraints)
logger.info("> Constraint equilibiration complete. Initialising next stage...")
if ions_present:
    _ = bs.nptrun_mem(timestep=2, no_steps=20000000, dump_time=10000, ref_temp=temp, ref_press=press, continuation="yes", constrain=False, lipid=lipid)
else:
    _ = bs.nptrun_mem(timestep=2, no_steps=20000000, dump_time=10000, ref_temp=temp, ref_press=press, continuation="yes", constrain=False, ions=False, lipid=lipid)

subprocess.call("%s %s grompp -f npt_noconstrain.mdp -c npt_constrain.gro -p topol.top -o npt_noconstrain.tpr"%(gmx_cmd, quiet), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -s npt_noconstrain.tpr"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -s npt_noconstrain.tpr"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv confout.gro npt_noconstrain.gro", shell=True)


# Production
logger.info("> Equilibiration complete. Beginning Production...")
if ions_present:
    _ = bs.nptrun_mem(timestep=timestep, no_steps=no_steps, dump_time=dump_time, ref_temp=temp, ref_press=press, continuation="yes", constrain=False, equib=False, lipid=lipid)
else:
    _ = bs.nptrun_mem(timestep=timestep, no_steps=no_steps, dump_time=dump_time, ref_temp=temp, ref_press=press, continuation="yes", constrain=False, equib=False, ions=False, lipid=lipid)

subprocess.call("%s %s grompp -f npt_STID.mdp -c npt_noconstrain.gro -p topol.top -o npt_STID.tpr"%(gmx_cmd, quiet), shell=True)

if n_proc != 1:
    subprocess.call("mpirun -np %i gmx_mpi %s mdrun -ntomp %i %s -s npt_STID.tpr"%(n_proc, quiet, ntomp, gpu_cmd), shell=True)
else:
    subprocess.call("%s %s mdrun -s npt_STID.tpr"%(gmx_cmd, quiet), shell=True)

subprocess.call("mv confout.gro npt_STID.gro", shell=True)
subprocess.call("mv traj.trr npt_STID.trr", shell=True)

# Parsing phase
logger.info("> Gromacs run successfully completed, beginning parsing the data")
logger.info("> Converting trajectory to multipdb")
subprocess.call("echo '1' | %s %s trjconv -f npt_STID.trr -s npt_STID.tpr -skip %i -o parse_traj.trr -pbc mol -ur compact"%(gmx_cmd, quiet, skip), shell=True)

logger.info("> Cleaning up, and moving the transmembrane region back to the origin...")
subprocess.call(current_p + "/parse_mem.sh", shell=True)

#name_tmp = pdb.split('.')[0]
#subprocess.call("mv sim.pdb %s_sim.pdb"%(name_tmp), shell=True)
logger.info("> Gromacs run completed successfully. Check %s_sim.pdb for the final result and input to the next stage, note it is probably a good idea to change the name of this file."%(name_tmp))
