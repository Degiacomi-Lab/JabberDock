# -*- coding: utf-8 -*-
# Copyright (c) 2014 EPFL (Ecole Polytechnique federale de Lausanne)
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
# Author : Giorgio Tamo, giorgio.tamo@epfl.ch
# Web site : http://lbm.epfl.ch


import numpy as np
from Protein import Protein
from copy import deepcopy
import time

#in CG segmemts, 0=receptor, 1=ligand

class AssemblyHeteroMultimer:
    def __init__(self,list_of_structures):
        #self.ligand_file = ligand_file
        #self.receptor_file = receptor_file
        #self.cg_atoms=cg_atoms

        self.list_of_structures = list_of_structures[0]
        self.list_of_names = list_of_structures[1]


    ###################
    #PRIVATE METHODS###
    ###################

    def _move_mobile_structure_to_origin(self, index):
        #get the center of geometry
        xyzCenter = np.mean(self.structure_list_coords[index],axis=0)
        self.structure_list_coords[index] -= xyzCenter
        #if len(self.cg_atoms)>0:
            #xyzCenter_cg=np.mean(self.cg_atoms[self.cg_atoms[:,5]==1,2:5],axis=0)
            #self.cg_atoms[self.cg_atoms[:,5]==1,2:5]-= xyzCenter_cg

    def _translate(self,index,coords):
        self.structure_list_coords[index] += np.array([coords[0],coords[1],coords[2]])

        #if len(self.cg_atoms)>0:
            #self.cg_atoms[self.cg_atoms[:,5]==1,2:5]+= np.array([self.coords[0],self.coords[1],self.coords[2]])

    def _rotation(self, index, coords):
    #angle in numpy need to be given in rad -> rad = deg * pi/180
        alpha = np.radians(coords[3])
        beta = np.radians(coords[4])
        gamma = np.radians(coords[5])
    #rotation around x
    #|1     0                0         |
    #|0     np.cos(alpha)      -np.sin(alpha)|
    #|0     np.sin(alpha)   np.cos(alpha) |
        Rx = np.array([[1,0,0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
        rotation = np.dot(Rx,np.dot(Ry,Rz))
        #multiply rotation matrice with each atom of the monomer
        self.structure_list_coords[index] = np.dot(self.structure_list_coords[index],rotation)
        #if len(self.cg_atoms)>0:
            #self.cg_atoms[self.cg_atoms[:,5]==1,2:5]=np.dot(self.cg_atoms[self.cg_atoms[:,5]==1,2:5],rotation)
            ##rotate cg ligand in cg matrix


    ################
    #PUBLIC METHODS#
    ################

    def place_all_mobile_structures (self, pos):
        self.structure_list_coords = []
        coords_array = deepcopy(pos[: (len(self.list_of_structures)-1)*6 ]) # only take the new coordinates of the rigid monomers #  len -1 
        coords_array.shape = ((len(self.list_of_structures)-1), 6) # len -1
        coords = []
        for i in xrange(0,len(self.list_of_structures),1):
            self.structure_list_coords.append(deepcopy(self.list_of_structures[i].monomer.data[:,5:8]))
             # here will have to select the first 6 when flexibility will be included

            # NB. the first structure is the biggest and is fixed
            # so you place all the others
            if i > 0:
                coords = coords_array[i-1] # i-1 because the first structure is fixed therefore has no new positions to be assigned
                self._move_mobile_structure_to_origin(i)
                self._rotation(i, coords)
                self._translate(i, coords)

    #def place_ligand(self, coords):

        #self.coords = coords
        #self.ligand = []
        #self.ligand = deepcopy(self.ligand_file.data[:,5:8])
        #self.receptor = []
        #self.receptor = deepcopy(self.receptor_file.data[:,5:8])

        ####print "start: %s"%self.cg_atoms[90,2:5]
        #self._move_ligand_to_origin()
        #self._rotation()
        #self._translate()

    def atomselect_structure(self,structure_index, chain,resid,atom,get_index=False):
        [m,index]=self.list_of_structures[structure_index].monomer.atomselect(chain,resid,atom,True)
        atoms = self.structure_list_coords[structure_index][index]

        if get_index==True:
            return [atoms, index]
        else:
            return atoms

    def atomselect (self, name_of_structure, chain,resid,atom):
#               print name_of_structure
        atom = self.atomselect_structure(self.list_of_names[name_of_structure], chain,resid,atom)
        return atom
    #def atomselect_receptor(self,chain,resid,atom,get_index=False):
    #[m,index]=self.receptor_file.atomselect(chain,resid,atom,True)
    #atoms = self.receptor[index]
    #if get_index==True:
        #return [atoms, index]
    #else:
        #return atoms

    #def get_width(self):
    #print ">> before to get the width here is the multimer = %s"%(self.multimer)
    #    maxXYZ = self._get_max_from_multimer()
    #    minXYZ = self._get_min_from_multimer()
    #print "get width took %s"%(end-start)
    #Simple way to calculate but take too much time
    #self.BigAtomArray = np.reshape(self.multimer,(-1,3))
    #self.MaxXYZ = np.amax(self.BigAtomArray,axis=0)
    #self.MinXYZ = np.amin(self.BigAtomArray,axis=0)

    #    return(maxXYZ[0]-minXYZ[0])

    #def get_height(self):
    #    maxXYZ = self._get_max_from_multimer()
    #    minXYZ = self._get_min_from_multimer()
    #    return(maxXYZ[2]-minXYZ[2])

    def get_structure_xyz(self, index):
        return self.structure_list_coords[index]

    #def get_receptor_xyz(self):
    #return self.receptor

    def distance(self,atom1,atom2):
        atom1np = np.array(atom1[0])
        atom2np = np.array(atom2[0])
        diff = atom1np - atom2np
        return np.sqrt(np.dot(diff,diff))

    def write_PDB(self,outname):

        f_out=open(outname,"w")

#               self.ligand_file.set_xyz(self.ligand)

        data_list = []

        chain_converter = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z')


        for structure_index in xrange(0,len(self.list_of_structures),1):
            self.list_of_structures[structure_index].monomer.set_xyz(self.structure_list_coords[structure_index])
            data_list.append(self.list_of_structures[structure_index].monomer.mapping(self.list_of_structures[structure_index].monomer.data))

        #map intergers to characters from ligand data
#               data_list=self.ligand_file.mapping(self.ligand_file.data)

        for j in xrange(0,len(data_list),1):

            for i in xrange(0,len(data_list[j]),1):
                for name, number in self.list_of_names.iteritems():
                    if number == j:
                        break

                #create and write PDB line # name is replacing chain converter here
                l=(data_list[j][i][0],data_list[j][i][1],data_list[j][i][2], name[0] ,data_list[j][i][4],data_list[j][i][5],data_list[j][i][6],data_list[j][i][7],data_list[j][i][8],data_list[j][i][9],data_list[j][i][10])
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)

            f_out.write("TER\n")


        #map intergers to characters from receptor data
#               data_list=self.receptor_file.mapping(self.receptor_file.data)
#
#               for i in xrange(0,len(data_list),1):
#                       #create and write PDB line
#                       l=(data_list[i][0],data_list[i][1],data_list[i][2],"R",data_list[i][4],data_list[i][5],data_list[i][6],data_list[i][7],data_list[i][8],data_list[i][9],data_list[i][10])
#                       L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
#                       f_out.write(L)


        f_out.close()

# ------------- FUNCTION to export a pdb file so that can be converted to density map ----------------------

    def create_PDB_for_density_map (self, procNo):

        #TODO: maybe you re going to have to create a density map for each processor?
        
        # all the monomer have already have had their xyz set already in the fitness class
        f_out=open("simulated_map"+str(procNo)+".pdb","w")
        data_list = []

        for structure_index in xrange(0,len(self.list_of_structures),1):
            self.list_of_structures[structure_index].monomer.set_xyz(self.structure_list_coords[structure_index])
            data_list.append(self.list_of_structures[structure_index].monomer.mapping(self.list_of_structures[structure_index].monomer.data))

        #map intergers to characters from ligand data
#               data_list=self.ligand_file.mapping(self.ligand_file.data)

        for j in xrange(0,len(data_list),1):

            for i in xrange(0,len(data_list[j]),1):
                for name, number in self.list_of_names.iteritems():
                    if number == j:
                        break

                #create and write PDB line # name is replacing chain converter here
                l=(data_list[j][i][0],data_list[j][i][1],data_list[j][i][2], name[0] ,data_list[j][i][4],data_list[j][i][5],data_list[j][i][6],data_list[j][i][7],data_list[j][i][8],data_list[j][i][9],data_list[j][i][10])
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)

            #f_out.write("TER\n")
        f_out.close()

# ------------- FUNCTIONS for the Coarse grain forcefield implementation -----------------------------------

    def get_CG_coords(self):
        return self.cg_atoms


    def get_CG_ligand(self):
        ###print "return: %s"%self.cg_atoms[90,2:5]
        return  self.cg_atoms[self.cg_atoms[:,5]==1]


    def get_CG_receptor(self):
        return  self.cg_atoms[self.cg_atoms[:,5]!=1]
