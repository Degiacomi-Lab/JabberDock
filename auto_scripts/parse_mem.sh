#!/bin/bash

loc=$(dirname $(echo "$0"))

# Script to parse out just the protein molecule of interest (from the membrane / water system)
# then move simulation by the transmembrane region to the origin (i.e. COM of TM to origin)
vmd -dispdev text -e ${loc}/align_lipid.tcl
python ${loc}/shift_tm.py
