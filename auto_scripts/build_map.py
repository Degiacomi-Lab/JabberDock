#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess
import numpy as np
import biobox as bb

# Current path destination:
current_p = str(os.path.dirname(os.path.realpath(__file__)))

# Create arguments based on user input / default 
parser = argparse.ArgumentParser(description='build map input parameters')
parser.add_argument('-i', metavar="input", required=True, help='Input simulation multipdb (e.g. sim.pdb from gmx_run)')
parser.add_argument('-ff_dat', metavar="forcefield", required=True, help='Forcefield charge data used to run the MD. This is a required option, but one is avaliable in biobox/classes/amber03 to correspond to the amber03 forcefield used in gmx_run.')
parser.add_argument('-o', metavar="outname", default='map', required=False, help='Name of prefix density / dipole file (the latter only if requested - default is your receptorname_map)')
parser.add_argument('-ts', metavar="time_start", default=15, required=False, help='Frame from which to begin building the STID maps from (default is 15)')
parser.add_argument('-te', metavar="time_end", default=105, required=False, help='Frame from which to stop building the STID maps from (default is 105)')
parser.add_argument('-t', metavar="temp", default=310.15, required=False, help='Temperature simulation was run at (default is 310.15K)')
parser.add_argument('-p', metavar="press", default=1.0, required=False, help='Pressure simulation was run at (default is 1.0 bar)')
parser.add_argument('-vox', metavar="vox_in_window", default=3., required=False, help='Number of voxels to include from the surrounding space into the central voxel (default is 3)')
parser.add_argument('-w', metavar="width", default=1., required=False, help='Width of a voxel (default is 1 Ang., units in ang.)')
parser.add_argument('-dm', metavar="dip_map", default=1, required=False, help='Whether you want the dipole map printing too (required for postprocessing - default is true). It needs an input of either 0 or 1 (false or true)')
parser.add_argument('-ac', metavar="ac", default=1, required=False, help='You are converting from an amber forcefield. Needs 0 or 1 - False or True (Default = True)')
parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
parser.add_argument('-eps', metavar="eps_ext", default=54., required=False, help='External Dielectric used in the calculation of the internal, default is 54.')

args = vars(parser.parse_args())

pdb = str(args["i"])
ff = str(args["ff_dat"])
outname = str(args["o"])
time_start = int(args["ts"])
time_end = int(args["te"])
temp = float(args["t"])
press = float(args["p"]) * 1e+5 # convert to Pa
vox = float(args["vox"])
width = float(args["w"])
dm = int(args["dm"])
eps_ext = float(args["eps"])
ac = int(args["ac"]) 

# log initial command
cmd_output_string = sys.argv
cmd_output = ""
for i in range(len(cmd_output_string)):
    cmd_output = cmd_output + cmd_output_string[i] + " "

if outname == 'map':
    protein = pdb.split('.')[0]
    outname = protein + '_' + 'map'

# create logger
logname = "build_map_run.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("build_map")
fh = logging.FileHandler(logname)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)

if not args["v"]:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

logger.info(cmd_output)

if dm < 0 or dm > 1:
    raise Exception("ERROR: Use a 1 or a 0 to indicate if you want a dipole map or not")
elif dm == 1:
    dip_map = True
else:
    dip_map = False

# Run error checks if time is too short / long
if time_start >= time_end:
    raise Exception('ERROR: time_end must be at least 1 frame more than time_start')

if time_end != -1:
    time_end += 1 # shift by 1 to account for python counting

# Begin building grid map
logger.info("> Building the voxel map...")

M = bb.Molecule()
M.import_pdb(pdb)

window_size = vox * width # resolution, in angstroms
crds = M.coordinates[time_start:time_end, :]

if ac == 1:
    M_pqr = M.pdb2pqr(ff=ff)
else:
    M_pqr = M.pdb2pqr(ff=ff, amber_convert=False)

M.assign_atomtype()
mass = M.get_mass_by_atom()

# Find the max / min boundaries of our box and add buffer regions to use desired resolutions and get a complete picture
maxmin = []
diff = []
remain = []
buff = []
buffmaxmin = []
for i in range(3):
    # 3 for the three cartisian axes, we want do define our coordinate system and reference frame as the FIRST in the sim.
    maxmin.append(np.amin(crds[:,:,i]))
    maxmin.append(np.amax(crds[:,:,i]))
    diff.append(maxmin[2*i+1] - maxmin[2*i])
    remain.append(diff[i] % window_size)

    buff.append(remain[i] / 2.)
    buffmaxmin.append(maxmin[2*i]-buff[i] - 2* window_size )
    buffmaxmin.append(maxmin[2*i+1]+buff[i] + 2* window_size )

# From minimum of buffer to maximum of buffer + voxel shift - res. () In other words, the start site of every coordinate.
# This is for _item in the loop above, where we analyse what atoms are present between this and the next (+ res)
x_range = np.arange(buffmaxmin[0], buffmaxmin[1]+width-window_size, width)   
y_range = np.arange(buffmaxmin[2], buffmaxmin[3]+width-window_size, width)  
z_range = np.arange(buffmaxmin[4], buffmaxmin[5]+width-window_size, width)

spec_vol = window_size**3 * 10**(-27) # -27 for docking # Convert nm^3 to metres^3
orig = np.array([x_range, y_range, z_range]) + window_size / 2.
a = np.array((len(x_range), len(y_range), len(z_range)))
min_crds = np.ndarray.tolist(np.asarray((buffmaxmin[0] + buff[0], buffmaxmin[2] + buff[1], buffmaxmin[4] + buff[2])) + window_size / 2)

# Build the dipole map
logger.info("> Grid built, creating a dipole map")
dipoles  = M.get_dipole_map(orig = orig, pqr = M_pqr, time_start = time_start, time_end = time_end, write_dipole_map = dip_map, vox_in_window = vox, resolution = width, fname = outname + ".tcl")

# Build the STID map
logger.info("> Dipole map written, building the STID map...")
M.get_dipole_density(dipole_map = dipoles, orig = orig, min_val = min_crds, vox_in_window = vox, V = spec_vol, outname = outname + ".dx", T = temp, P = press, resolution = width, epsilonE=eps_ext)

logger.info("> Creating a single PDB structure centered in your maps")

M.write_pdb(outname + '.pdb', conformations=[int((time_start+time_end)/2)])

logger.info("> STID map building complete. Check %s.pdb, %s.dx and %s.tcl for output"%(outname, outname, outname))
