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


import numpy as np
from Protein import Protein
from copy import deepcopy
import time

class Multimer:
    def __init__(self,pdb):
        #store a Protein structure
        self.pdb = pdb


    def _move_extreme_point_to_origin(self):
        #get the extreme point on the x axis and move the atom corresponding to it to the origin
        xyzMaxIndex = np.argmax(self.monomer,axis=0)
        maxAtom = self.monomer[xyzMaxIndex[0]]
        self.monomer -= np.array([maxAtom[0],maxAtom[1],0])


    def _move_monomer_to_origin(self):
        xyzCenter = np.mean(self.monomer,axis=0)
        self.monomer -= xyzCenter


    def _translate(self):
        self.monomer -= np.array([self.radius,0,0])


    def _rotation(self):
        #angle in numpy need to be given in rad -> rad = deg * pi/180
        alpha = np.radians(self.pos[0])
        beta = np.radians(self.pos[1])
        gamma = np.radians(self.pos[2])
        #rotation autour axe x
        #|1     0                0         |
        #|0     np.cos(alpha)      -np.sin(alpha)|
        #|0     np.sin(alpha)   np.cos(alpha) |
        Rx = np.array([[1,0,0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
        rotation = np.dot(Rx,np.dot(Ry,Rz))
        #multiply rotation matrix with each atom of the monomer
        self.monomer = np.dot(self.monomer,rotation)


    def _circular_symmetry(self):
        #create the number of monomer required for the multimer and rotate them around the z axis
        for i in xrange(0,self.degree,1):
                
            #number of degree to rotate
            self.angle = np.radians(i*(360/self.degree))
            #print "angle %s for multimer %s"%(angle,i)
            self.Rz = np.array([[np.cos(self.angle), -(np.sin(self.angle)), 0], [(np.sin(self.angle)), (np.cos(self.angle)), 0], [0,0,1]])
            self.multimer.append(np.dot(self.monomer,self.Rz))
            

    def _get_max_from_multimer(self):
        arrayMax = []
        #calculate the position of each monomer maximum and store them in an array
        for m in xrange(0,self.degree,1):
            arrayMax.append(np.amax(self.multimer[m],axis=0))

        #get the maximum among all monomer
        maxXYZ = np.amax(arrayMax,axis=0);
        return maxXYZ


    def _get_min_from_multimer(self):
        arrayMin = []
        #calculate the position of each monomer minimum and store them in an array
        for m in xrange(0,self.degree,1):
            arrayMin.append(np.amin(self.multimer[m],axis=0))

        #get the minimum among all monomer
        minXYZ = np.amin(arrayMin,axis=0);
        return minXYZ


    def create_multimer(self, degree, radius, pos):

        self.degree = degree
        self.radius = radius
        self.pos = pos
        self.multimer = []
        self.monomer = deepcopy(self.pdb.data[:,5:8])

        self._move_monomer_to_origin()
        self._rotation()
        self._move_extreme_point_to_origin()
        self._translate()
        self._circular_symmetry()


    def multimer_to_origin(self):
        m=self.get_multimer_xyz()
        center = np.mean(m,axis=0)
        for i in xrange(0,self.degree,1):
            self.multimer[i][:,0] -= center[0]
            self.multimer[i][:,1] -= center[1]
            self.multimer[i][:,2] -= center[2]


    def z_to_origin(self):
        zCenter = np.mean(self.multimer[0],axis=0)[2]
        for i in xrange(0,self.degree,1):
            self.multimer[i][:,2] -= zCenter


    def z_shift(self,z):
        for i in xrange(0,self.degree,1):
            self.multimer[i][:,2] += z


    def z_rotation(self,angle):
        #rotation of the whole assembly around around z
        theta = np.radians(angle)
        Rz = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0,0,1]])
        for i in xrange(0,self.degree,1):
            M=self.multimer[i]
            self.multimer[i] = np.dot(M,Rz)


    def atomselect(self,unit,chain,resid,atom,get_index=False):
        [m,index]=self.pdb.atomselect(chain,resid,atom,True)
        atoms = self.multimer[unit-1][index]

        if get_index==True:
            return [atoms, index]
        else:
            return atoms


    def get_width(self):
        maxXYZ = self._get_max_from_multimer()
        minXYZ = self._get_min_from_multimer()

        return(maxXYZ[0]-minXYZ[0])


    def get_height(self):
        maxXYZ = self._get_max_from_multimer()
        minXYZ = self._get_min_from_multimer()
        return(maxXYZ[2]-minXYZ[2])


    def get_multimer_uxyz(self):
        return self.multimer


    def get_multimer_xyz(self):
        multimerxyz = []
        for i in xrange(0,len(self.multimer),1):
            multimerxyz.extend(self.multimer[i])
        return multimerxyz


    def distance(self,atom1,atom2):
        #atom1np = np.array(atom1[0])
        #atom2np = np.array(atom2[0])
        #diff = atom1np - atom2np
        #return np.sqrt(np.dot(diff,diff))

        d=[]
        for i in xrange(0,len(atom1),1):
            d.append(np.sqrt(np.sum((atom2-atom1[i])**2,axis=1)))
        
        dist=np.array(d)
        return np.min(d)



    def write_PDB(self,outname):

        chain_converter = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z')

        f_out=open(outname,"w")

        for j in xrange(0,self.degree,1):
            self.pdb.set_xyz(self.multimer[j])

            #map intergers to characters from input data (default: all the protein)
            data_list=self.pdb.mapping(self.pdb.data)

            for i in xrange(0,len(data_list),1):
                #create and write PDB line
                l=(data_list[i][0],data_list[i][1],data_list[i][2],chain_converter[j],data_list[i][4],data_list[i][5],data_list[i][6],data_list[i][7],data_list[i][8],data_list[i][9],data_list[i][10])
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)

            f_out.write("TER\n")

        f_out.close()








