#!/bin/bash

loc=$(dirname $(echo "$0"))

# hack topol issue introducted by combining slipid with amber
# Sometimes gromacs is... picky and will use A2 instead of A etc. You will get an error:
# ERROR 29258 [file topol_Protein_chain_A.itp, line 146783]:....
# etc if this is the case. Simply change A2 to A or whatever in this script (corresponding to whats now in the folder you're running in), to resolve this
# e.g. change B to B2 if you have more than a monomer transmembrane protein (naming convention changes)
python ${loc}/hack_topol.py topol_Other_chain_B.itp 5 2
python ${loc}/hack_topol.py topol_Other_chain_A2.itp 5 2
mv topol_Other_chain_B_edit.itp topol_Other_chain_B.itp
mv topol_Other_chain_A2_edit.itp topol_Other_chain_A2.itp
