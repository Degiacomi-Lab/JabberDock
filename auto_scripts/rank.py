#!/usr/bin/env python

import argparse
import sys, os
from copy import deepcopy
import logging
import subprocess
import numpy as np
from scipy.stats import rankdata

############################################
################ Parameters ################
############################################

# Create parameters based on user arguments

# Need number of samples, input style

parser = argparse.ArgumentParser(description='Parse input parameters')
parser.add_argument('-d', metavar='data_file', required=False, default='models/additional_files/model_solutions.dat', help="Data file containing the scores. This can be either the raw scores output from dock or from dipole_rerank")
parser.add_argument('-o', metavar='output_file', required=False, default='models/ranked_scores.dat', help="Output re-ranked scores. Note that this will contain an additional column to your input, the original line number, as these correspond to the model numbers for the output pdb files")
parser.add_argument('-v', action="store_true", help='Verbose (I want updates!)')
args = vars(parser.parse_args())

in_file = str(args['d'])
out_file = str(args['o'])

#test to see if files exist
if not os.path.isfile(in_file):
    raise Exception("ERROR: File %s not found!"%in_file)

############################################

# create logger
logname = "rank.log"
if os.path.isfile(logname):
    os.remove(logname)

logger = logging.getLogger("rank")
fh = logging.FileHandler(logname)
ch = logging.StreamHandler()
logger.addHandler(fh)
logger.addHandler(ch)

if not args["v"]:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

# log initial command
cmd_output_string = sys.argv
cmd_output = ""
for i in range(len(cmd_output_string)):
    cmd_output = cmd_output + cmd_output_string[i] + " "

logger.info(cmd_output)

logger.info("> Loading in %s"%(in_file))
data = np.loadtxt(in_file)

dat_shape = data.ndim
numbers = np.arange(0, np.shape(data)[0], 1, dtype=int)

if dat_shape == 1:
    num_data = np.concatenate((data[:, None], numbers[:, None]), axis=1)
elif dat_shape == 2:
    num_data = np.concatenate((data, numbers[:, None]), axis=1)
else:
    raise Exception("ERROR: Your input data file has a strange file structure. It needs a specific number of columns to be accepted")

# Ranking function
def rank_scores(data, index):
    """
    Take data and desired column to order by, and return a corresponding ranked list
    :param data: Input data in the form of a 1D array
    :param index: Column index to use to rank the data
    """
    
    # Get column
    col = data[:, index]

    # Get list of ranks
    ind = np.argpartition(col, -len(col))[-len(col):]
    list_data = data[ind]
    col_list = col[ind]

    # Get ranked list
    rank = np.absolute(rankdata(col_list, method='dense') - len(col_list))

    # Reorder data based on ranks
    ranked_data = np.asarray([x for _, x in sorted(zip(rank, list_data))])

    return ranked_data


if dat_shape == 1:
    logger.info("> Picked up the dipole re-ranked file. Reranking accordingly...")
    ranked_data = rank_scores(data = num_data, index = 0)

elif dat_shape == 2:
    logger.info("> Picked up the dock output file. Reranking accordingly...")
    ranked_data = rank_scores(data = num_data, index = 7)
    ranked_data = ranked_data[:, [7,-1]]

else:
    raise Exception("ERROR: File is not dipole re-ranked or dock output file. It has %i columns when it should have 1 or 9"%(dat_shape))

logger.info("> Ranking complete, outputting ranked data to %s. The numbers in the last column indicate the corresponding model numbers of the complex PDB file output by dock."%(out_file))

f = open(out_file, "w")
f.write("SCORE  MODEL_NO")
for line in ranked_data:
    f.write("%.5f  %i\n"%(line[0], int(line[1])))

f.close()
