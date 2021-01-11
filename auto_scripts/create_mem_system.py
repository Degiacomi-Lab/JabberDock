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

"""Build the membrane system"""

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments
parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-i', metavar="pdb", required=True, help='Input PDB file')
parser.add_argument('-ff', metavar="ff", default='amber14sb_slipid', required=False, help='Forcefield reference name (default amber03_slipid. Note, MUST be merged forcefield with lipid parameters present)')
parser.add_argument('-rep', metavar="repair", default="T", required=False, help='Repair the structure using modeller? (default T, True). N.B., requires Modeller and biopython module installed')
args = vars(parser.parse_args())

pdb = str(args["i"])
ff = str(args["ff"])
repair = str(args["rep"])

if "amber" in ff.lower() and "slipid" in ff.lower():
    amber_con = "amber"
else:
    amber_con = "na"

if repair == "T" or repair == "True":
    repair_struc = "repair"
else:
    repair_struc = "na"

subprocess.call(current_p + "/create_embed.sh %s %s %s"%(pdb, amber_con, repair_struc), shell=True)

print(">> System setup complete, refer to gmx_mem_run.py to continue the process")
