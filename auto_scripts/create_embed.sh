#!/bin/bash

# setup files and run preliminary gromacs routines (and fix forcefield errors)
loc=$(dirname $(echo "$0"))

input_file=${1}
ff=${2} # using amber with slipid? (otherwise following hacks won't work)
rep=${3}

# begin by fixing any missing atoms in the residues with modeller
if [ ${rep} == "repair" ]; then
    bash ${loc}/repair_chain.sh $(cut -d'.' -f1 <<<${input_file})
fi

# use packmol-memgen on desired prealigned protein (N.B this adds no salt)
packmol-memgen -l POPE:POPE -r 1:1 --nocounter -p ${input_file} --preoriented --nloop 80 --tolerance 2.4 --random --nloop_all 120 

# correct for POPE resname and atomnames
if [ ${ff} == "amber" ]; then
    echo "> Hacking resnames to accomodate amber and slipid usage..."
    bash ${loc}/replace_lipid.sh bilayer_${input_file} fixed.pdb
else
    mv bilayer_${input_file} fixed.pdb
fi

# use vmd to remove silly added protein hydrogens
vmd -dispdev text -e ${loc}/rem_prot_h.tcl

# correct fixed_trim.pdb by adding a TER between protein and lipid (vmd removes TER)
if [ ${ff} == "amber" ]; then
    sed -i '/\bC32 POPEA   1\b/ i TER' fixed_trim.pdb
fi

# get resid of TM region saved for later shifting of coordinates
python ${loc}/trans_res.py 
