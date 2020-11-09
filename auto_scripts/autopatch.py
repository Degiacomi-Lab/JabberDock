#! /usr/bin/python

# given a PDB file missing some residues and its expected FASTA sequence, create a complete model using MODELLER
# if import of modeller fails, first source modellerX/bin/modpy.sh, e.g.:
# source /home/xdzl45/bin/modeller9.19/bin/modpy.sh

import sys
import os
import re
from textwrap import wrap
from modeller import *
from modeller.automodel import * 

class ChainError():
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

def _prepare_alignment(fbasename):
    
    fasta_fname = fbasename+".fasta"
    
    f = open(fasta_fname, "r")
    f1 = f.readlines()
    f.close()
    P1_pos = []
    for i in range(len(f1)):
        if (">" in f1[i]):
            P1_pos.append(i)
    seq_name = f1[P1_pos[0]].split('|')[0].replace('>', '')
    
    w = open("alignment.seg","w+")
    
    if len(P1_pos) == 1:
        # we have a homodimer
        seq = f1[1][:-1] + '/' + f1[1][:-1] + '*\n'
    elif len(P1_pos) == 2:
        # sequence written style needs the '\', we also want to remove the \n
        seq = f1[1][:-1] + '/' + f1[3][:-1] + '*\n'
    else:
        raise ChainError('More than two chains detected!')
    
    # always 7 lines to write to
    w.write('>P1;%s\n'%(fbasename.lower()))
    w.write('sequence::::::::-1.00:-1.00\n')
    w.write(seq)
    w.write('\n')
    w.write('>P1;%s\n'%(fbasename))
    w.write('structureX:%s:FIRST:@:END :@:::-1.00:-1.00\n'%(fbasename))
    w.write('*')
    w.close()

def autopatch(fbasename, gap_cutoff=8):

    #pdb_out = "%s_PATCHED.pdb"%fbasename; the output pdb file name (if successful, empty otherwise) 
    pdb_out = ""
    seq_name = fbasename+'.seq'
    try:

        _prepare_alignment(fbasename)
        _pdb2seq(fbasename)
        _fasta2pir(fbasename)
        _full_multi_align(fbasename)
        _trim_align("multi_align.ali")
        patch_status = _gap_check("trimmed_align.ali", gap_cutoff)

        if patch_status == "yes":
            pdb_out = _patch_model(fbasename, seq_name)
        else:
            pdb_out = ""

    except Exception as e:
        print("ERROR: %s"%e)
        return ""

    myfiles = ['%s.seq'%fbasename, '%s.pir'%fbasename, 'alignment.seg', 'multi_align.ali',
               'alignment.seg.ali', 'trimmed_align.ali', 'family.mat',
               '*.ini', '*.rsr', '*.sch', '%s.V*'%fbasename.lower(), '%s.D*'%fbasename.lower()]
    
    for m in myfiles:
        try:
            os.system("rm %s &> /dev/null"%m)
        except:
            continue    

    return pdb_out

###################################################################

#step 1a. pir format of AA from pdb
def _pdb2seq(fbasename):
    env = environ()
    mdl = model(env, file=fbasename)
    aln = alignment(env)
    aln.append_model(mdl, align_codes=fbasename)
    aln.write(file=fbasename+'.seq')

#step 1b. pir from complete AA fasta 
def _fasta2pir(fbasename):
    env = environ()
    a = alignment(env, file=fbasename+".fasta", alignment_format='FASTA')
    a.write(file=fbasename+'.pir', alignment_format='PIR')

#for multi chain
def _full_multi_align(fbasename): # align the seq from *.pdb and RCSB_fasta

    log.verbose()
    env = environ()
    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(env,
              # file with template codes and target sequence
              alnfile  = 'alignment.seg',
              # PDB codes of the templates (structureX section) in alignment.seg
              knowns   = fbasename,
              # code of the target (sequence section) in alignment.seg
              sequence = fbasename.lower())
    a.auto_align()                      # get an automatic alignment

    # replace '-' with '/' after chain in the alignment.seg.ali file.
    align_file = "alignment.seg.ali"
    f=open(align_file, "r")
    f1 = f.readlines()
    f.close()
    P1_pos = []

    #identify positions of different sequances (P1) blocks
    for i in range(len(f1)):
        if ("P1;" in f1[i]):
            P1_pos.append(i)
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]

    #remove empty lines and split into 2 blocks: seq and structure
    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1+2:sec_2]])
    AA_seq = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_2+2:]])

    chain_break_seq = AA_seq.index('/')
    AA_struc_new = AA_struc[:chain_break_seq] + '/' + AA_struc[chain_break_seq + 1:]

    AA_struc_new = wrap(AA_struc_new, 75) #split after 75 characters 
    AA_seq_new = wrap(AA_seq, 75)

    AA_struc_new = '\n'.join([str(elem) for elem in AA_struc_new])
    AA_seq_new = '\n'.join([str(elem) for elem in AA_seq_new])
    
    struc_sec = f1[sec_1] + f1[sec_1+1] + AA_struc_new + "\n"
    seq_sec = f1[sec_2] + f1[sec_2+1] + AA_seq_new + "\n"
    f = open("multi_align.ali", 'w')
    f.writelines(struc_sec)
    f.writelines(seq_sec)
    f.close()

#step 2. add sequence name to 2nd line; copy the pir contents and structure info into alignment.seg; align sequences and generate model
def _full_align(fbasename):
    pir_fname = fbasename+'.pir'
    seq_fname = fbasename+'.seq'
    f = open(pir_fname, "r")
    f1 = f.readlines()
    f.close()
    #f1 = [x.rstrip() for x in f1]
    for i in range(len(f1)):
        if ("P1;" in f1[i]):
            seq_name = f1[i][4:8]
        if ("sequence:" in f1[i]):
            seq_pos = i 
    P1 = ">P1;"+seq_name+"\n"
    seq_line = "sequence:"+seq_name+":::::::-1.00:-1.00\n"
    AA_block = ''.join([str(elem) for elem in f1[seq_pos+1:]])
    seq_block = P1 + seq_line + AA_block
    f = open(pir_fname, 'w')
    f.writelines(seq_block)
    f.close()

    myCmd_A = 'cat %s %s > alignment.seg'%(pir_fname, seq_fname)
    os.system(myCmd_A)

    env = environ()
    env.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(env,
                  # file with template codes and target sequence
                  alnfile  = 'alignment.seg',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name) 
    a.auto_align()   # get an automatic alignment (alignment.seg.ali)
    return seq_name

#step 3. trim the alignment by removing gaps for missing residues at the termini of the structure
def _trim_align(align_file):
    align_file = "alignment.seg.ali"
    f=open(align_file, "r")
    f1 = f.readlines()
    P1_pos = []

    #identify positions of different sequances (P1) blocks
    for i in range(len(f1)):
        if ("P1;" in f1[i]):
            P1_pos.append(i)
    f.close()
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]


    #remove empty lines and split into 2 blocks: seq and structure
    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1+2:sec_2]])

    start_gaps = re.findall('^[-]+', AA_struc)
    end_gaps = re.findall('[-]+$', AA_struc)

    if len(start_gaps)>0:
        start_gaps_len = len(start_gaps[0])
    else:
        start_gaps_len = 0

    if len(end_gaps)>0:
        end_gaps_len = len(end_gaps[0])
    else:
        end_gaps_len = 0

    AA_struc_len = len(AA_struc)
    AA_struc_new = AA_struc[start_gaps_len : AA_struc_len - end_gaps_len]+"*"
    AA_struc_new = wrap(AA_struc_new, 75) #split after 75 characters 
    
    AA_seq = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_2+2:]])
    AA_seq_len = len(AA_seq)
    AA_seq_new = AA_seq[start_gaps_len : AA_seq_len - end_gaps_len]+"*"
    AA_seq_new = wrap(AA_seq_new, 75) 
    
    AA_struc_new = '\n'.join([str(elem) for elem in AA_struc_new])
    AA_seq_new = '\n'.join([str(elem) for elem in AA_seq_new])
    
    struc_sec = f1[sec_1] + f1[sec_1+1] + AA_struc_new + "\n"
    seq_sec = f1[sec_2] + f1[sec_2+1] + AA_seq_new + "\n"
    f = open("trimmed_align.ali", 'w')
    f.writelines(struc_sec)
    f.writelines(seq_sec)
    f.close()

#step 4. Check if any gap is more than cutoff length in the trimmed_align.ali and if so set patch_status = "no"
def _gap_check(align_file, gap_cutoff):
    patch_status = "yes"
    align_file = "trimmed_align.ali"
    f=open(align_file, "r")
    f1 = f.readlines()
    P1_pos = []
    for i in range(len(f1)):
        if ("P1;" in f1[i]):     #identify positions of different sequances (P1) blocks
            P1_pos.append(i)
    f.close()
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]

    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1    +2:sec_2]])
    #check for gap_lengths in AA_struc
    struc_gaps = re.findall('[-]+', AA_struc)
    for gaps in range(len(struc_gaps)):
        gap_len = len(struc_gaps[gaps])
        if gap_len > gap_cutoff:
            print("pdb not patched. Long gap: "+ str(gap_len))
            patch_status = "no"

    return patch_status

#step 5. build missing residues
def _patch_model(fbasename, seq_name):
    print(">> patching model...")
    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(env,
                  # file with template codes and target sequence
                  alnfile  = 'trimmed_align.ali',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = fbasename.lower(), 
                  assess_methods = (assess.DOPE, assess.GA341))     
    a.md_level = refine.fast #very_fast, fast, slow, very_slow, slow_large, refine
    #repeat whole cycle twice and do not stop unless obj. func > 1e6
    #a.repeat_optimization = 2
    a.max_molpdf = 1e6
    a.make()

    pdb_out = "%s_PATCHED.pdb"%fbasename
    myCmd_mv = 'mv %s.B99990001.pdb %s'%(fbasename.lower(), pdb_out)
    os.system(myCmd_mv)
    
    return pdb_out

#######################################################

if __name__ == "__main__":

    # fname should be a basename: expect to find both a .pdb and a .fasta file with that name
    fbasename = sys.argv[1]

    # gap should be an integer, if user fails to provide one, use a dega
    try:
        gap = int(sys.argv[2])
    except:
        print("setting up default gap size=10")
        gap = 10

    foutname = autopatch(fbasename, gap)
    if foutname == "":
        print("autopatch failed")
    else:
        print("saved patched file %s"%foutname)


