#! /usr/bin/python


import numpy as np
import time


def surface_complementarity(dens1, dens2, delta1, delta2, o1, o2):
    return -1

def chris_magic(dens1, dens2, delta1, delta2, o1, o2):
    return -1


#load all data
receptor = np.load("receptor.npy")
good_ligand = np.load("good_ligand.npy")
bad_ligand = np.load("bad_ligand.npy")

receptor_delta = np.load("receptor_delta.npy")
good_ligand_delta = np.load("good_ligand_delta.npy")
bad_ligand_delta = np.load("bad_ligand_delta.npy")

receptor_origin= np.array([])
good_ligand_origin= np.array([])
bad_ligand_origin= np.array([])


# call original surf. comp. evaluation function
start = time.time()
score1 = surface_complementarity(receptor, good_ligand, receptor_delta, good_ligand_delta, receptor_origin, good_ligand_origin)
score2 = surface_complementarity(receptor, bad_ligand, receptor_delta, bad_ligand_delta, receptor_origin, bad_ligand_origin)
end = time.time()
print("> Surface complemetarity:")
print("> time: %s"%(end-start))
print("> good score: %s"%score1)
print("> bad score: %s"%score2)


# call Chris' magical solution
start = time.time()
score1 = chris_magic(receptor, good_ligand, receptor_delta, good_ligand_delta, receptor_origin, good_ligand_origin)
score2 = chris_magic(receptor, bad_ligand, receptor_delta, bad_ligand_delta, receptor_origin, bad_ligand_origin)
end = time.time()
print("\n> Chris magic:")
print("> time: %s"%(end-start))
print("> good score: %s"%score1)
print("> bad score: %s"%score2)

