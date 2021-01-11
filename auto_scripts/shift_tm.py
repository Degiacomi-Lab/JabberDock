import biobox as bb
import numpy as np

# Shift a simulations coordinates by the previously recorded TM region's COM to the origin
res = np.load('TM_res.npy') # resid of TM region

A = bb.Molecule()
A.import_pdb('sim.pdb')

A_atoms = A.atomselect("*", res, "*")

COM = np.mean(A_atoms, axis=0)

A.translate(-COM[0], -COM[1], -COM[2])
A.guess_chain_split()
A.write_pdb('sim.pdb')
