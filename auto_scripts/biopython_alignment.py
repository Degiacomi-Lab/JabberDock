# run as python3
# seq alignment better than match_residue in biobox
import biobox as bb
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 

def seq_knowledge(res):
    # return the amino acid code for the input sequence
    return {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'PHE': 'F',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
    }.get(res, 'X')    # X is default if not found

def convert_res(res):
    # do loop to get correct sequence for input res
    seq = ""
    for r in res:
        true_res = r[-3:] # remove C or N prefix
        seq += seq_knowledge(true_res)
    
    return seq

def get_alignment(Rb, Ru):
    """
    get alignment of two PDBs, based on the biopython algorithm and two input PDBs, name1 and name2
    """
    
    # get CA for alignment
    Rb_idxs = Rb.atomselect("*", "*", "CA", get_index=True)[1]
    Ru_idxs = Ru.atomselect("*", "*", "CA", get_index=True)[1]
    
    Rb_sub = Rb.get_subset(Rb_idxs)
    Ru_sub = Ru.get_subset(Ru_idxs)
    
    # convert to sequence code
    Rb_seq = convert_res(Rb_sub.data["resname"])
    Ru_seq = convert_res(Ru_sub.data["resname"])
    
    R_alignment = pairwise2.align.globalxx(Rb_seq, Ru_seq)
    txt = format_alignment(*R_alignment[-1])
    
    return txt

def condense_sequence(seq, matches):
    """
    condence a seq down since the PDB doesn't know about the gaps in terms of the resid numbering.
    i.e. we want the index (or boolean) of the resids we need for the pdb, not the extended aligned seq
    """
    blanks = np.asarray(["-" in a for a in seq])
    #seq = np.asarray(list(seq))
    
    AA = np.invert(blanks)
    return matches[AA] # return boolean matches of correct length with removed blank spaces

def match_residue(name1, name2):
    """
    Match residues of two pdbs, similar to bb.match_residue, but better with biopython
    """

    # Import structures
    R1 = bb.Molecule()
    R2 = bb.Molecule()
    
    R1.import_pdb(name1) 
    R2.import_pdb(name2)

    # Assign atomtypes just in case
    R1.assign_atomtype()
    R2.assign_atomtype()

    # we need to align to the pdb with more residues (otherwise we might get funky matching answers)
    # check the length of each first
    R1_CA = R1.atomselect("*", "*", "CA", get_index=True)[1]
    R2_CA = R2.atomselect("*", "*", "CA", get_index=True)[1]

    if len(R1_CA) > len(R2_CA):
        # get alignment
        ali = get_alignment(R1, R2)
        R1_seq, lines, R2_seq = ali.split('\n')[:3]
    else:
        ali = get_alignment(R2, R1)
        R2_seq, lines, R1_seq = ali.split('\n')[:3]

    # now only get resid matches which works for both
    matches = np.asarray(["|" in a for a in lines])

    matches1 = condense_sequence(R1_seq, matches)
    matches2 = condense_sequence(R2_seq, matches)

    # need to match where there is a line and a letter in each
    return matches1, matches2



