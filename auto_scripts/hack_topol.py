"""
This script replaces the angle and dihedral functional types in an input file to the proposed arguments
e.g. if you want all angle types in topol_Other_chain_X.itp to be 5 and all dihedrals to be 4:
python topol_hack.py input.itp 5 4
N.B. this is inefficient and dirty (too many nested loops) - it could also be done in an awk script much more efficiently
"""

import sys, os

filename = str(sys.argv[1])
angles = int(sys.argv[2])
dihedrals = int(sys.argv[3])
outname = filename.split('.')[0] + '_edit.' + filename.split('.')[1]

lines = open(filename).read().splitlines()

with open(outname, 'w') as write_file:

    editmarker_angle = False
    editmarker_dihedral = False
    # loop through itp file and write all lines as normal until we come across angles or dihedrals
    for line in lines:
        if editmarker_angle:
            if line == '':
                write_file.write(line + '\n')
                editmarker_angle = False
            elif len(line) > 35:
                write_file.write(line + '\n') # write back in comment
            else:
                writeline = line[:-2] + '%i '%(angles) # all but last two characters (for function)
                write_file.write(writeline + '\n')
        elif editmarker_dihedral:
            if line == '':
                write_file.write(line + '\n')
                editmarker_dihedral = False
            elif len(line) > 87:   # catch for proper dihedrals (want to edit improper dihedrals only)
                write_file.write(line + '\n')
                editmarker_dihedral = False
            elif len(line) < 87 and len(line) > 32:
                write_file.write(line + '\n') # write back in comment
            else:
                writeline = line[:-2] + '%i '%(dihedrals)
                write_file.write(writeline + '\n')
        else:
            write_file.write(line + '\n')

        if line[:10] == '[ angles ]':
            editmarker_angle = True
        elif line[:14] == '[ dihedrals ]':
            editmarker_dihedral = True
        else:
            continue

write_file.close()

