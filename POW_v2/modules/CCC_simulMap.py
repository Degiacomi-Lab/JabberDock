# -*- coding: utf-8 -*-
import numpy as np
import copy
import scipy.integrate
from scipy import ndimage
import scipy.ndimage.filters
import subprocess

def make_simulated_map (pdb_file, rank, voxel_space, resolution ):
#    subprocess.call("./kercon "+str(pdb_file)+" simulated_map.sit "+str(voxel_space)+" "+str(resolution))
    subprocess.call("FilaSitus/bin/kercon "+str(pdb_file)+" simulated_map"+str(rank)+".sit "+str(voxel_space)+" "+str(resolution)+" > out", shell=True)


def compute_corr(exp_map, simul_map, resolution):

    matrix1file = open(exp_map)
    matrix2file = open(simul_map)
    matrix1 = []
    matrix2 = []
    tempomatrix1 = []
    tempomatrix2 = []
    header1 = []
    header2 = []
    header_flag1 = 1
    header_flag2 = 1

    # reading the .srt file
    while (1):
        line = matrix1file.readline()

        if (line and header_flag1):
            noNewLine = line.replace('\n','')
            header1.append(noNewLine.split(" "))
            header_flag1 = 0
    #        sequence += [space.split(" ")]
            #sequence[-1].remove('\n')
        elif (line and header_flag1 == 0 ):
            noNewLine = line.replace('\n','')
            tempomatrix1.append(noNewLine.split(' '))
            # extracting everything elements of the temporary matrix and cleaning:
            for strings in tempomatrix1:
                for element in strings:
                    if len(element) != 0:
                        matrix1.append(element)
            tempomatrix1 = []
        else:
            break

    tempoMat1 = np.float32(matrix1)

    while (1):
        line = matrix2file.readline()

        if (line and header_flag2):
            noNewLine = line.replace('\n','')
            header2.append(noNewLine.split(" "))
            header_flag2 = 0
    #        sequence += [space.split(" ")]
            #sequence[-1].remove('\n')
        elif (line and header_flag2 == 0 ):
            noNewLine = line.replace('\n','')
            tempomatrix2.append(noNewLine.split(' '))
            # extracting everything elements of the temporary matrix and cleaning:
            for string in tempomatrix2:
                for element in string:
                    if len(element) != 0: 
                        matrix2.append(element)
            tempomatrix2 = []
        else:
            break

    tempoMat2 = np.float32(matrix2)
    
    
    
    
    if resolution < 10:
        # leave the maps unchanged
        numpyMat = copy.deepcopy(tempoMat1)
        numpyMat2 = copy.deepcopy(tempoMat2)
    
    elif resolution >= 10:
        # apply laplacian filter on both density maps:
        numpyMat = scipy.ndimage.filters.laplace(tempoMat1)
        numpyMat2 = scipy.ndimage.filters.laplace(tempoMat2)
    
    
    # Situs based cross correlation coefficient:
    if len(numpyMat) <= len(numpyMat2):
        numerator = scipy.integrate.quad(lambda r: (numpyMat[r] * numpyMat2[r]), 0, len(numpyMat))
        Deno1 = scipy.integrate.quad(lambda r: (numpyMat[r] * numpyMat[r]), 0, len(numpyMat))
        Deno2 = scipy.integrate.quad(lambda r: (numpyMat2[r] * numpyMat2[r]), 0, len(numpyMat))
    else:
        numerator = scipy.integrate.quad(lambda r: (numpyMat[r] * numpyMat2[r]), 0, len(numpyMat2))
        Deno1 = scipy.integrate.quad(lambda r: (numpyMat[r] * numpyMat[r]), 0, len(numpyMat2))
        Deno2 = scipy.integrate.quad(lambda r: (numpyMat2[r] * numpyMat2[r]), 0, len(numpyMat2))
        
    if Deno2[0] != 0 and Deno1[0] != 0: # because denominator cannot be equal to 0
        corr = numerator[0] / ( (Deno1[0] ** (0.5)) * (Deno2[0] ** (0.5)) )
    else:
        corr = 0.0
    print ">>> CCC = "+str(corr)
    return corr
