import biobox as bb
import numpy as np

# get the residues in the TM region for later adjustment of COM of transmembrane to origin at end of simulation period
zmax = 19.4 # bounds of the TM region (as found through an average of 14 bilayers)
zmin = -19.4

A = bb.Molecule()
A.import_pdb('fixed_trim.pdb')

A_atoms, A_idxs = A.atomselect("*", "*", "CA", get_index=True)

A_atoms_z = A_atoms[:,2]

TM_loc = ~np.logical_or(A_atoms_z > zmax, A_atoms_z < zmin)
A_idxs_TM = A_idxs[TM_loc]

res = np.asarray(A.get_subset(A_idxs_TM).data['resid'])
np.save('TM_res', res)
