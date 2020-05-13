# -*- coding: utf-8 -*-
# Copyright (c) 2012 EPFL (Ecole Polytechnique federale de Lausanne)
# Laboratory for Biomolecular Modeling, School of Life Sciences
#
# POW is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# POW is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with POW ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Author : Matteo Degiacomi, matteothomas.degiacomi@epfl.ch
# Web site : http://lbm.epfl.ch


#! /usr/bin/env python


import numpy as np
import scipy
from scipy.cluster.vq import *

import sys, os
from copy import deepcopy
import subprocess
import re

from MDAnalysis import *
import MDAnalysis.analysis.align


# read in a topology and a trajectory and re,turn a matrix of atomic coordinates (also write a pdb of first frame)
def loader(trajectory, topology, align="yes",selection="protein"):

   #check input consistency
    if os.path.isfile(topology)!=True :
        if rank==0:
            print "ERROR: topology file not found!"
        sys.exit(1)

    if os.path.isfile(trajectory)!=True :
        if rank==0:
            print "ERROR: trajectory file not found!"
        sys.exit(1)

    print ">> loading trajectory..."
    #load trajectory and save reference structure
    universe=Universe(topology,trajectory)
    protein_sel = universe.selectAtoms(selection)
    protein_sel.write("protein.pdb")

    #align trajectory on reference frame
    if align=="yes":
        print ">> aligning trajectory on protein..."
        reference=Universe("protein.pdb")
        MDAnalysis.analysis.align.rms_fit_trj(universe,reference,select=selection, filename="aligned.dcd")
        universe=Universe(topology,"aligned.dcd")

    #extract protein's atoms coords
    print ">> extracting atomic coordinates..."
    protein = universe.trajectory.timeseries(universe.selectAtoms(selection),format="fac")
    all_coords=protein.reshape(len(protein),len(protein[0])*len(protein[0][0])).transpose()

    return all_coords


#run PCA analysis on C alpha of a trajectory and write main projections of proj_coordinates.dat file
def pca(trajectory, topology, ratio):

    #check input consistency
    if os.path.isfile(topology)!=True :
        if rank==0:
            print "ERROR: topology file not found!"
        sys.exit(1)

    if os.path.isfile(trajectory)!=True :
        if rank==0:
            print "ERROR: trajectory file not found!"
        sys.exit(1)

    if ratio>=1 or ratio<=0:
        if rank==0:
            print "ERROR: ratio of simulation energy should be a number in ]0;1[!"
        sys.exit(1)

    #extract protein's atoms coords
    print ">> preparing data for PCA analysis..."
    universe=Universe(topology,trajectory)
    ca_sel = universe.selectAtoms("name CA")
    ca = universe.trajectory.timeseries(ca_sel,format="fac")
    coords=ca.reshape(len(ca),len(ca[0])*len(ca[0][0])).transpose()
    print ">>> found %s frames and %s degrees of freedom"%(len(coords[0]),len(coords[:,0]))

    #check whether conditions for PCA analysis are good
    if len(coords[0])<len(coords[:,0]):
        print "ERROR: trajectory length should be greater than the system's degrees of freedom (3N)!"
        sys.exit(1)

    #compute displacement matrix (removing mean pos from atom pos in coords matrix)
    disp=deepcopy(coords.astype(float))
    for i in xrange(0,len(disp),1):
        disp[i]-=np.mean(disp[i])

    #compute covariance matrix, eigenvalues and eigenvectors
    print ">> computing covariance matrix..."
    covariance=np.cov(disp)
    print ">> extracting eigenvalues and eigenvectors (might take few minutes)..."
    [eigenval,eigenvec]=np.linalg.eig(covariance)

    #compute representative number of eigenvectors according to desired ratio (user provided)
    print "\n   nb., eigenvalue, cumulative ratio"
    print "   ---------------------------------"
    cumulative=0
    cnt=0
    for i in xrange(0,len(eigenval),1):
        cumulative+=eigenval[i]
        print "   %s, %s, %s"%(i+1, eigenval[i], cumulative/np.sum(eigenval))
        if(cumulative/np.sum(eigenval)>ratio):
            cnt=i+1
            break

    #compute projection of trajectory on the significant number of eigenvectors
    #lines = n-th eigenvector component, columns = simulation frame
    print "\n>> projecting trajectory on %s eigenvectors..."%cnt
    p=[]
    for i in xrange(0,cnt,1):
        p.append(np.dot(eigenvec[:,i],coords))

    proj=np.array(p)

    #output trajectory projection into file
    f = open("proj_coordinates.dat", 'w')
    for i in xrange(0,len(proj.transpose()),1):
        ln=str(proj.transpose()[i])+"\n"
        f.write(re.sub("\]","",re.sub("\[","",ln)))
    f.close()
    return proj,eigenvec


#projections clustering space with kmeans algo
# c=number of required clusters (future: decide dynamically this value)
def clusterize(proj):

    c=1
    print ">> clustering frames in the projections space..."
    res,idx=kmeans(proj.transpose(),c)
    code,min_dist=vq(proj.transpose(),res)
    frames_out=[]
    a=np.ma.array(min_dist)
    for x in xrange(0,c,1):
        a.mask=(code!=x)
        frames_out.append(a.argmin())
    print ">>> frame %s selected as protein configurations centroid"%frames_out[0]

    return frames_out[0]


#deform a structure along a linear combination of eigenvectors
def deform_structure(X_deform,deform_coeffs,eigenvec):
    for n in xrange(0,len(eigenvec),1):
        X_deform+=deform_coeffs[n]*eigenvec[:,n]
    return X_deform


#extract from coordinates the one being the closest to the centroid
def get_nearest(all_coords, proj, centroid):
    #find which frame in trajectory is the closest to the desired deformation
    #the first centroid is kept as reference (should replace 0 by an iterator, when using multiple centroids)
    pos=proj[:,centroid]+deform_coeffs
    code,min_dist=vq(proj.transpose(),np.array([pos]))
    target_frame=min_dist.argmin()

    return all_coords[:,target_frame]
