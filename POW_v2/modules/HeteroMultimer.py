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


from Default import Parser as R
from Default import Space as S
from Default import Postprocess as PP

import numpy as np
import os, sys
from copy import deepcopy
from scipy.cluster.vq import *

from Protein import Protein
import AssemblyHeteroMultimer as A
import flexibility_new as F
import Multimer as M

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

import CCC_simulMap as CCC

import ClusterAndDraw as CnD
import wx

class Parser(R):

    def __init__(self):

        self.add('rmsdSelect','rmsd_select','array str',"NA")
        self.add('constraint','constraint','str',"NA")
        self.add('energy','energy_type','str',"vdw")
        self.add('detectClash','detect_clash','str',"on")
        self.add('target','target','array float',np.array([]))
        self.add('mode','mode','str',"seed")
        self.add('mixingWeight','mix_weight','float', 0.2)

        # style of the assembly, rigid or flexible
        self.add('assembly_style','assembly_style','str',"rigid")

        # flexibility flags
        self.add('topology', 'topology', 'array str', 'NA')
        self.add('trajectory','trajectory', 'array str','NA')
        self.add('trajSelection','trajselection','str','NA')
        self.add('ratio','ratio','float',0.9)
        self.add('align', 'align', 'str', 'no')
        self.add('projection','proj_file','str',"NA")

        # monomer flags
        self.add('monomer','monomer_file_name','array str', "NA")

        # post processing flags
        self.add('cluster_threshold','cluster_threshold','float',"NA")
        self.add('output_folder','output_folder','str',"result")

        # density map docking flag
        self.add('density_map_docking','map_dock_OnOff', 'str', 'no')
        self.add('density_map','density_map', 'str', 'NA')


    def check_variables(self):

        if self.cluster_threshold<0:
            print "ERROR: clustering threshlod should be greater than 0!"
            sys.exit(1)

        # GIORGIO_CODE check existence of pdb files for the monomers in the folder XXXXXXXXXXXXXXXXXXXXXX
        for nickname,pdb_file in self.monomer_file_name:
            if pdb_file != "NA" and nickname != "NA": #  and self.monomer_style=="rigid"
                tmp=os.path.abspath(pdb_file)
                if os.path.isfile(pdb_file)!=1 :
                    print "ERROR: monomer pdb file %s not found"%pdb_file
                    sys.exit(1)

        # CHECKING THE TRAJECTORIES AND THE TOPOLOGY FILE
        if self.assembly_style=="flexible":
            for nickname,top_file in self.topology:
                if top_file != "NA" and nickname != "NA": #  and self.monomer_style=="rigid"
                    tmp=os.path.abspath(pdb_file)
                    if os.path.isfile(top_file)!=1 :
                        print "ERROR: monomer pdb file %s not found"%top_file
                        sys.exit(1)

            for nickname,traj_file in self.trajectory:
                if traj_file != "NA" and nickname != "NA": #  and self.monomer_style=="rigid"
                    tmp=os.path.abspath(traj_file)
                    if os.path.isfile(traj_file)!=1 :
                        print "ERROR: monomer pdb file %s not found"%traj_file
                        sys.exit(1)

            # checking that every trajectory has a topology:
            if len(self.trajectory) != len(self.topology):
                print "ERROR: unequal number of topologies and trajectories"

        #check if number of target measures are provided
        if len(self.target) == 0 :
            print 'ERROR: target measures not specified!'
            sys.exit(1)

        #test if constraint file exists, and function constaint_check is available
        #try:
        exec 'import %s as constraint'%(self.constraint.split('.')[0])
        #except:
        #    print "ERROR: load of user defined constraint function failed!"
        #    sys.exit(1)
        try:
            constraint.constraint_check
        except AttributeError:
            print 'ERROR: constraint_check function not found in file %s'%self.constraint
            sys.exit(1)

        # check whether density flag is on or off
        if self.map_dock_OnOff == "on" :
            if self.density_map != "NA":
                self.map_dock_OnOff = True
            else:
                print 'ERROR: if density map flag is on, a density map should be present in the directory'
                sys.exit(1)
            print ">> Electron Density Map Docking Mode"
        elif self.map_dock_OnOff == "off":
            self.map_dock_OnOff = False
        else:
            print 'ERROR: density_map_docking should be either "on" or "off"'
            sys.exit(1)



class Structure():

    def __init__ (self):
            # initialising the values
        self.monomer = "NA" # the monomeric unit
        self.pdb_file_name = "NA"
        self.index_CA_monomer = "NA"
        self.flexibility = "NA"
        self.init_coords = "NA"

    def read_pdb (self, pdb):
        self.pdb_file_name = pdb
        self.monomer = Protein()
        self.monomer.import_pdb(pdb)
        self.init_coords = self.monomer.get_xyz()

    def compute_PCA (self, topology,trajectory,align,ratio,mode, proj_file):
        self.flexibility = F.Flexibility_PCA()
        self.flexibility.compute_eigenvectors(topology,trajectory,align,ratio,mode, proj_file)


    def setCoords (self):
        self.init_coords = self.monomer.get_xyz()


    #def get_index



class Data:

    index_ligand=[]
    index_receptor=[]
    cg_atoms=[]

    def __init__(self,params):

        self.structure_hash = {}
        self.structure_list = []
        volume_structure_hash = {} # this is used to get the biggest structure


        # GIORGIO_CODE create structure instances of the rigid monomers
        for nickName,pdb_file in params.monomer_file_name:
            # create instance of the structure class
            s = Structure()
            s.read_pdb (pdb_file)

            volume_structure_hash[len(s.monomer.get_xyz())] = [s, nickName]


        # create structure instance for the flexible monomers
        if params.assembly_style=="flexible":
            print ">> flexible docking requested for structures, launching PCA..."

            for nickName, traj_file in params.trajectory:
                try:
                    # get the topology file:
                    for nickName2, top_file in params.topology:
                        if nickName2 == nickName:
                            break

                    # create the structure and compute the PCA
                    s = Structure()
                    s.compute_PCA(top_file, traj_file, params.align, params.ratio, params.mode, params.proj_file)
                    s.read_pdb("protein.pdb")

                    volume_structure_hash[len(s.monomer.get_xyz())] = [s, nickName]

                except ImportError, e:
                    sys.exit(1)

                # TODO: work on the deform mode, but ask Matteo before
                if params.mode=="deform":
                    self.structure_ligand=Protein()
                    self.structure_ligand.import_pdb("protein.pdb")
                    self.ligand.import_pdb("CA.pdb")

        # getting the biggest structure and putting at the beginning so that it is fixed
        sorted_volumes = volume_structure_hash.keys()
        sorted_volumes.sort()
        sorted_volumes.reverse()

        for i in sorted_volumes:

            # insert the elements in a list
            self.structure_list.append( volume_structure_hash[i][0] ) # insert the structure
            self.structure_hash[volume_structure_hash[i][1]] = self.structure_list.index(volume_structure_hash[i][0])

        self.structure_list_and_name = [self.structure_list, self.structure_hash]
        print self.structure_list_and_name

        #LIGAND STRUCTURE
        #self.ligand = Protein()

#               if params.assembly_style=="flexible":
#                       print ">> flexible docking requested for ligand, launching PCA..."
#                       try:
#                               self.flex_ligand=F.Flexibility_PCA()
#                               self.flex_ligand.compute_eigenvectors(params.ligand_topology,params.ligand_trajectory,params.ligand_align,params.ligand_ratio,params.mode,params.ligand_proj_file)
#                               self.ligand.import_pdb("protein.pdb") # importing the middle structure
#                       except ImportError, e:
#                               sys.exit(1)
#
#                       if params.mode=="deform":
#                               self.structure_ligand=Protein()
#                               self.structure_ligand.import_pdb("protein.pdb")
#                               self.ligand.import_pdb("CA.pdb")

        #else:
            #load monomeric structure (the pdb file)
            #self.ligand.import_pdb(params.ligand_file_name)

        if params.energy_type=="vdw":
            self.CA_index_of_structures = self.get_index(["CA"])
            #[self.index_ligand,self.index_receptor]=self.get_index(["CA","CB"])

        # if the density map docking is on load the structure into data:
        if params.map_dock_OnOff:
            self.density_map_fileName = params.density_map
            
    def get_index(self,atoms=["CA","CB"]):

        #generate a dummy assembly and extract the indexes where atoms of interest are located
        # first create the numpy array containing all null translation and rotation for each of the mobile structures
        null_coordinate_array = np.zeros((len(self.structure_list)-1)*6)
        assembly = A.AssemblyHeteroMultimer(self.structure_list_and_name)
        assembly.place_all_mobile_structures(null_coordinate_array)

        #ligand_index=[]
        #receptor_index=[]
        index_of_all_structures = [] # this is going to be an array of arrays

        for aname in atoms:

            for structure_number in xrange(0,len(self.structure_list),1):
                index_of_all_structures.append([])

                #append indexes of an element in atoms list for all structures
                [m,index]=assembly.atomselect_structure(structure_number , "*","*",aname,True)
                for i in index:
                    index_of_all_structures[structure_number].append(i)

                ##append indexes of an element in atoms list for receptor
                #[m,index]=assembly.atomselect_receptor("*","*",aname,True)
                #for i in index:
                    #receptor_index.append(i)

        return index_of_all_structures



class Space(S):
    def __init__(self,params,data):

        len_flexi=0
        if params.assembly_style=="flexible":
            for structure in data.structure_list:
                if structure.flexibility != "NA":
                    len_flexi += len(structure.flexibility.eigenspace_size)


        len_rec=0
        #if params.receptor_style=="flexible":
            #len_rec=len(data.flex_receptor.eigenspace_size)

        len_rigid_dim = 6*(len(data.structure_list)-1)

        # for hetero-multimer assembly, given that every MOBILE (so exept the first one) protein has 6 degrees of freedom
        self.low=np.zeros(len_rigid_dim +len_flexi)
        self.high=np.zeros(len_rigid_dim +len_flexi)
        self.cell_size=np.zeros(len_rigid_dim +len_flexi)
        self.boundary_type=np.zeros(len_rigid_dim +len_flexi)


        #box size as given by all the structures dimensions
        first_min=np.min(data.structure_list[0].monomer.get_xyz(),axis=0)
        first_max=np.max(data.structure_list[0].monomer.get_xyz(),axis=0)

        distance_array = []

        for x in xrange (1,len(data.structure_list),1):
            distance_array.append(np.max(data.structure_list[x].monomer.get_xyz(),axis=0) - np.min(data.structure_list[x].monomer.get_xyz(),axis=0))

        distance_array = np.array(distance_array)

        summed_distances = np.sum(distance_array, axis = 0)

        box_min=first_min-(summed_distances)
        box_max=first_max+(summed_distances)


        if len(params.high_input)!=len(params.low_input):
            print "ERROR: boundaryMin and boundaryMax should have the same length!"
            sys.exit(1)

        #assign low boundaries
        if params.low_input!="NA" :

            if len(params.low_input)== len_rigid_dim :
                for i in xrange(0,len(params.low_input),1):
                    self.low[i]=params.low_input[i]
            else:
                print "ERROR: boundaryMin should contain 6 values (3 rotations, 3 translations)"
                sys.exit(1)
        else:
            print "WARNING: boundaryMin undefined, using default values"
            i = 0
            for x in xrange(0, len_rigid_dim ,1):
                if i < 3:
                    self.low[x] = box_min[i]
                    i += 1
                elif (i > 2) and (i != 6):
                    self.low[x] = 0.0
                    i+=1
                if i == 6:
                    i = 0


        #assign high boundaries
        if params.high_input!="NA" :

            if len(params.high_input)== len_rigid_dim:
                for i in xrange(0,len(params.high_input),1):
                    self.high[i]=params.high_input[i]
            else:
                print "ERROR: boundaryMax should contain 6 values (3 rotation, 3 translation)"
                sys.exit(1)
        else:
            print "WARNING: boundaryMax undefined, using default values"
            i = 0
            for x in xrange(0, len_rigid_dim ,1):
                if i < 3:
                    self.high[x] = box_max[i]
                    i += 1
                elif (i > 2) and (i != 6):
                    self.high[x] = 360.0
                    i+=1
                if i == 6:
                    i = 0

        # add all the flexible structures eigenvector fluctuations in the search space
        if params.assembly_style=="flexible":
            i = 0
            for structure in data.structure_list:
                if structure.flexibility != "NA":
                    for x in xrange(0, len(structure.flexibility.eigenspace_size),1):
                        self.low[len_rigid_dim+i]=-structure.flexibility.eigenspace_size[x]
                        self.high[len_rigid_dim+i]= structure.flexibility.eigenspace_size[x]
                        i += 1


        #add ligand eigenvector fluctuations in search space
#               for i in xrange(0,len_flexi,1):
#                       self.low[len_rigid_dim+i]=-data.flex_ligand.eigenspace_size[i]
#                       self.high[len_rigid_dim+i]=data.flex_ligand.eigenspace_size[i]



        #check boundary conditions consistency
        if len(self.low) != len(self.high):
            print 'ERROR: dimensions of min and max boundary conditions are not the same'
            sys.exit(1)
        if (self.low>self.high).any():
            print 'ERROR: a lower boundary condition is greated than a higher one'
            sys.exit(1)
        #define cell size
        self.cell_size=self.high-self.low

        #set boundary type (periodic for angles, repulsive for translation)
        if params.boundary_type=="NA":
            for i in xrange(0,len(self.low),1):
                self.boundary_type[i]=0
        elif params.boundary_type!="NA" and len(params.boundary_type)!=len(self.low):
            print 'ERROR: boundaries type inconsistent with system dimensions'
            print 'ERROR: %s dimensions given, but %s needed!'%(len(params.boundary_type),len(self.low))
            sys.exit(1)
        else:
            for i in xrange(0,len(self.low),1):
                self.boundary_type[i]=params.boundary_type[i]



class Fitness:
    def __init__(self,data,params):

        self.mode=params.mode
        self.map_docking_flag = params.map_dock_OnOff # do so because you want to pass this var to the function evaluate below
        
        # loading the reference/experimental density map file if flag is on
        if self.map_docking_flag:
            self.density_map_fileName = params.density_map
        
        #check if target exists
        try: params.target
        except NameError:
            print 'ERROR: target measures not found'
            sys.exit(1)
        self.target=params.target

        self.constraint=params.constraint.split('.')[0] # Xx constraint [rigidRandom].py
        #test if constraint file exists, and function constaint_check is available
        try:
            exec 'import %s as constraint'%(self.constraint)
        except ImportError, e:
            print "ERROR: load of user defined constraint function failed!"
            sys.exit(1)

        try: constraint.constraint_check
        except NameError:
            print 'ERROR: constraint_check function not found'

        #data to manipulate
        self.data=data

        self.assembly_style=params.assembly_style
        #self.receptor_style=params.receptor_style

        self.len_lig=0
        #if params.ligand_style=="flexible":
            #self.len_lig=len(self.data.flex_ligand.eigenspace_size)

        self.len_rec=0
        #if params.receptor_style=="flexible":
            #self.len_rec=len(self.data.flex_receptor.eigenspace_size)
        self.c1=params.mix_weight

    def evaluate(self,num,pos):

        exec 'import %s as constraint'%(self.constraint)
        import AssemblyHeteroMultimer as A

#               if ligand is flexible, select the most appropriate frame
        for structure in self.data.structure_list:
            if self.assembly_style=="flexible" and structure.flexibility != "NA":
                len_rigid_dim = 6*(len(self.data.structure_list)-1)
                i = 0

                deform_coeffs = self.coordinateArray[n][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]

                if self.mode=="seed":
                    pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                    code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                    target_frame=min_dist.argmin()
                    coords=structure.flexibility.all_coords[:,target_frame]
                    coords_reshaped=coords.reshape(len(coords)/3,3)
                    structure.monomer.set_xyz(coords_reshaped)
                else:
                    coords=structure.monomer.get_xyz()
                    coords_reshaped=coords.reshape(len(coords)*3)

                    for n in xrange(0,len(deform_coeffs),1):
                        coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                    structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                i += len(structure.flexibility.eigenspace_size)
            else:
                structure.monomer.set_xyz(structure.init_coords)


        # after getting the positions from PSO, create a new assembly according to those positions
        self.assembly = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
        self.assembly.place_all_mobile_structures(pos)
        
        # ------------------------------- DEFAULT FITNESS FUNCTION --------------------------
        
        if self.map_docking_flag == False:

            #if needed, compute error with respect of target measures
            distance=0
            if len(self.target)!=0:
    
                measure = constraint.constraint_check(self.data, self.assembly) # returns the distances between some ligand and receptor atoms
    
                if len(measure) != len(self.target) :
                    print 'ERROR: measure = %s'%measure
                    print 'ERROR: target measure = %s'%self.target
                    print 'ERROR: constraint file produced %s measures, but %s target measures are provided!'%(len(measure),len(self.target))
                    sys.exit(1)
    
                diff=self.target-np.array(measure)
                distance=np.sqrt(np.dot(diff,diff))
                #c1=0.1
    
            #compute system energy
            energy=0
            if len(self.data.CA_index_of_structures[0])>0:
    
                c1=self.c1 # was 0.2
                energy=self.interface_vdw()
                return c1*energy+(1-c1)*distance
                #return energy/len(self.data.index_ligand)+distance
            else:
                print "WHAT THE...???"
            #else:
            #       c1=0.001
            #       energy=self.measure_cg_energy(self.assembly,num)
            #       #fitness = coulomb+vdw+distance
            #       return c1*(energy[1]+energy[2])+(1-c1)*distance
            
        # --------------------------------- DENSITY MAP DOCKING ----------------------------
        
        elif self.map_docking_flag == True:
            
            # the coefficient for CCC, bigger one will give a heavier weight for density map docking
            c2 = 20 # try make vary from 1 to 100
            
            #proceed to calculate check for the geometry constraints
            distance=0
            if len(self.target)!=0:
    
                measure = constraint.constraint_check(self.data, self.assembly) # returns the distances between some ligand and receptor atoms
    
                if len(measure) != len(self.target) :
                    print 'ERROR: measure = %s'%measure
                    print 'ERROR: target measure = %s'%self.target
                    print 'ERROR: constraint file produced %s measures, but %s target measures are provided!'%(len(measure),len(self.target))
                    sys.exit(1)
    
                diff=self.target-np.array(measure)
                distance=np.sqrt(np.dot(diff,diff))
            
            # compute the systems energy
            energy=0
            if len(self.data.CA_index_of_structures[0])>0:
                c1=self.c1 
                energy=self.interface_vdw()

            fitness_score = c1*energy+(1-c1)*distance + c2 # the + 1 at the end is used to have a better score when good map dock fitting
            
            # -------------------- DENSITY MAP DOCKING FITNESS
            if fitness_score < (c2 + 10):
                
                resol = 15
                
                print ">>> Density map refinement rank "+str(rank)
                # create the pbd file to be transformed into the density map
                self.assembly.create_PDB_for_density_map(rank)
                #create the simulated density map
                CCC.make_simulated_map ("simulated_map"+str(rank)+".pdb", rank, 1, resol )
                #compare the two density maps and extract their cross correlation coefficient:
                ccc =  CCC.compute_corr(self.density_map_fileName, "simulated_map"+str(rank)+".sit", resol)
                
                #return the score of final function:
                return c1*energy+(1-c1)*distance+ c2*(1 - ccc)
            
            else:
                return fitness_score


    def measure_target(self):
        
        #measure constraints
        measure = constraint.constraint_check(self.assembly)

        if len(measure) != len(self.target) :
            print 'ERROR: measure = %s'%measure
            print 'ERROR: target measure = %s'%self.target
            print 'ERROR: constraint file produced %s measures, but %s target measures are provided!'%(len(measure),len(self.target))
            sys.exit(1)

        #measure distance within target values and obtained values
        diff=self.target-np.array(measure)
        distance=np.sqrt(np.dot(diff,diff))
        return distance


    def interface_vdw(self):

        epsilon=1.0
        sigma=4.7
        cutoff=12.0
        energy=0

        # for Heteromultimer assembly, you need to compute the energy of every structure against each other:
        for structure1_index in xrange (0,len(self.data.structure_list), 1):
            for structure2_index in xrange (structure1_index,len(self.data.structure_list), 1):
                d=[]
                if structure1_index == structure2_index:
                    pass

                else:
                    m1=self.assembly.get_structure_xyz(structure1_index)[self.data.CA_index_of_structures[structure1_index]]
                    m2=self.assembly.get_structure_xyz(structure2_index)[self.data.CA_index_of_structures[structure2_index]]

                    #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms
                    #m1=self.assembly.get_ligand_xyz()[self.data.index_ligand]
                    #m2=self.assembly.get_receptor_xyz()[self.data.index_receptor]

                    #extract distances of every atom from all the others
                    for i in xrange(0,len(m1),1):
                        d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

                    dist=np.array(d)

                    #detect interfacing atoms (atom couples at less than a certain cutoff distance
                    couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues


                    for i in xrange(0,len(couples[0]),1):
                        d=dist[couples[0,i],couples[1,i]]
                        energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
                        
        return energy



class Postprocess(PP):

    def __init__(self,data,params):
        self.data=data
        self.params=params


        self.len_lig=0
        #if params.assembly_style=="flexible":
            #self.len_lig=len(self.data.flex_ligand.eigenspace_size)

        self.len_rec=0
        #if params.receptor_style=="flexible":
            #self.len_rec=len(self.data.flex_receptor.eigenspace_size)

        #load constraint file
        self.constraint=params.constraint.split('.')[0]
        try:
            exec 'import %s as constraint'%(self.constraint)
        except ImportError, e:
            print "ERROR: load of user defined constraint function failed!"
            sys.exit(1)

        try:
            constraint.constraint_check
        except NameError:
            print 'ERROR: constraint_check function not found'


    #clustering according to rmsd of solutions in search space

    def run(self):
        if rank == 0:

            #create output directory for generated PDB
            self.OUTPUT_DIRECTORY=self.params.output_folder
            if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
                os.mkdir(self.OUTPUT_DIRECTORY)

            #use superclass method to filter acceptable solutions
            self.log=self.select_solutions(self.params) # -> the result is in fact the self.filter_log already
            print ">> %s solutions filtered"%len(self.log[:,1])
            if len(self.log[:,1])==0:
                return

            self.coordinateArray = deepcopy(self.log)   #[:, 0:len(self.log[0,:])].astype(float)
            self.dummyMatrix = np.empty(len(self.coordinateArray)**2)
            self.dummyMatrix.fill(100)
            self.distanceMatrix = self.dummyMatrix.reshape(len(self.coordinateArray),len(self.coordinateArray))
            self.dummyMatrix = []


            # variables to sliece the matrix into equal portions
            total_size = (len(self.coordinateArray)**2)/2
            binNo = size
            indexes_per_bin = total_size / binNo

            soustractor = 1
            indexBinHash = {}
            accumulator = 0

            rankIterator = 0
            lowBoundary = 0


            # getting the sliced indexes
            for i in xrange(0, len(self.distanceMatrix),1):
                array_len = len(self.distanceMatrix[i]) -  soustractor
                accumulator += array_len

                if accumulator > indexes_per_bin:
                    indexBinHash[rankIterator] = [lowBoundary, i]

                    # change the parameters
                    rankIterator += 1
                    lowBoundary = i
                    # empty the accumulator
                    accumulator = 0


                soustractor += 1

            if lowBoundary < i:
                indexBinHash[binNo-1] = [lowBoundary, i]


            print ">> Starting distance matrix creation:"
            print ">> Clustering best solutions..."

        else:
            self.distanceMatrix = None
            self.coordinateArray = None
            indexBinHash = None

        #synchronize all processers
        comm.Barrier()
        self.distanceMatrix=comm.bcast(self.distanceMatrix,root=0)
        self.coordinateArray=comm.bcast(self.coordinateArray,root=0)
        indexBinHash=comm.bcast(indexBinHash,root=0)
        comm.Barrier()


        exec 'import %s as constraint'%(self.constraint)



        #clusters_file=open("%s/dist_matrix.dat"%self.params.output_folder,"w") # Xx this where you write the solution file

        #generate a dummy multimer and extract the indexes of C alpha
#               null_coordinate_array = np.zeros((len(self.data.structure_list)-1)*6)
#               assembly = A.AssemblyHeteroMultimer(self.data.structure_list)
#               assembly.place_all_mobile_structures(null_coordinate_array)

        # Get the CA indexes already computed from the data class
        self.CA_index_of_all_structures = self.data.CA_index_of_structures


        #[m,index]=assembly.atomselect_of_structures(1,"*","*","CA",True) # -> extracting indexes of CA

        #load the monomeric structure positions (needed for resetting atom position after displacement)
#               s = Protein()
#               s.import_pdb(self.params.pdb_file_name)
#               coords=s.get_xyz()

        if len(self.coordinateArray) > (size *3):

            #----------------------------- first create the rmsd matrix
            # creating variables to check for status of clustering of process 0
            if rank == 0:
                repetitions = indexBinHash[rank][1] - indexBinHash[rank][0]
                totalIterations = len(self.coordinateArray) * repetitions
                counter = 0
                printresent = 1 # those are used not to repeat the state of the clustering
                printPast = 0

            counter = 0

            #synchronize all processes (get current timestep and repeat from swarm state)
            pieceOfCoordinateArray = np.array([])


            if rank in indexBinHash.keys():
            #Starting the creation with 2 loops
                for n in xrange(indexBinHash[rank][0],len(self.coordinateArray),1):

                    if n == indexBinHash[rank][1]:
                        break
                    for m in xrange (n,len(self.coordinateArray),1):
                        # make sure you are not using the same structures against themselves
                        if n == m:
        #                                       # add a "wrong" distance in the matrix to only have half the matrix
                            pass
                        else:

                            # --------------------------------- MODIFY THE FLEXIBLE STRUCTURES FOR THE 1ST ASSEMBLY AND SET COORDS
                            for structure in self.data.structure_list:
                                if self.params.assembly_style=="flexible" and structure.flexibility != "NA":
                                    len_rigid_dim = 6*(len(self.data.structure_list)-1)
                                    i = 0

                                    deform_coeffs = self.coordinateArray[n][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]

                                    if self.params.mode=="seed":
                                        pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                                        code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                                        target_frame=min_dist.argmin()
                                        coords=structure.flexibility.all_coords[:,target_frame]
                                        coords_reshaped=coords.reshape(len(coords)/3,3)
                                        structure.monomer.set_xyz(coords_reshaped)
                                    else:
                                        coords=structure.monomer.get_xyz()
                                        coords_reshaped=coords.reshape(len(coords)*3)

                                        for n in xrange(0,len(deform_coeffs),1):
                                            coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                                        structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                                    i += len(structure.flexibility.eigenspace_size)
                                else:
                                    structure.monomer.set_xyz(structure.init_coords)



                            # ------------------- CREATING 1ST ASSEMBLY
                            assembly1 = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
                            assembly1.place_all_mobile_structures(self.coordinateArray[n][:len(self.coordinateArray[n])-1])
                            # get the coordinates of all the structures to get the coordinates of Assembly
                            coordinate_of_assembly1_structures = []
                            for structure_index in xrange(0,len(self.data.structure_list), 1):
                                coordinate_of_assembly1_structures.append(assembly1.get_structure_xyz(structure_index))
                            m1 = np.concatenate((coordinate_of_assembly1_structures),axis=0)

                            # --------------------------------- MODIFY THE FLEXIBLE STRUCTURES FOR THE 2ND ASSEMBLY
                            for structure in self.data.structure_list:
                                if self.params.assembly_style=="flexible" and structure.flexibility != "NA":
                                    len_rigid_dim = 6*(len(self.data.structure_list)-1)
                                    i = 0

                                    deform_coeffs = self.coordinateArray[m][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]


                                    if self.params.mode=="seed":
                                        pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                                        code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                                        target_frame=min_dist.argmin()
                                        coords=structure.flexibility.all_coords[:,target_frame]
                                        coords_reshaped=coords.reshape(len(coords)/3,3)
                                        structure.monomer.set_xyz(coords_reshaped)
                                    else:
                                        coords=structure.monomer.get_xyz()
                                        coords_reshaped=coords.reshape(len(coords)*3)

                                        for n in xrange(0,len(deform_coeffs),1):
                                            coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                                        structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                                    i += len(structure.flexibility.eigenspace_size)
                                else:
                                    structure.monomer.set_xyz(structure.init_coords)

                            # ------------------- CREATING 2ND ASSEMBLY
                            assembly2 = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
                            assembly2.place_all_mobile_structures(self.coordinateArray[m][:len(self.coordinateArray[m])-1])
                            # get the coordinates of all the structures to get the coordinates of Assembly
                            coordinate_of_assembly2_structures = []
                            for structure_index in xrange(0,len(self.data.structure_list), 1):
                                coordinate_of_assembly2_structures.append(assembly2.get_structure_xyz(structure_index))
                            m2 = np.concatenate((coordinate_of_assembly2_structures),axis=0)



                            # calculate RMSD between the 2
                            rmsd=self.align(m1,m2) # --> comes from Default.Postprocess.align()
                            self.distanceMatrix[n][m] = rmsd

                            if rank == 0:
                                counter += 1.0
                                printPresent = int((counter / totalIterations) * 100)
                                if (printPresent%10) == 0 and printPresent != printPast:
                                    print "> ~"+str( printPresent )+" % structures clustered "
                                    printPast = printPresent

                pieceOfCoordinateArray = self.distanceMatrix[indexBinHash[rank][0]:indexBinHash[rank][1],:]
#                print " Clustering process "+str(rank)+" finished"

            comm.Barrier()
            pieces = comm.gather(pieceOfCoordinateArray,root=0)
            comm.Barrier()

            if rank == 0:
                self.distanceMatrix = []
                for elem in pieces:
                    if len(elem) < 2:
                        pass
                    else:
                        for arrays in elem:
                            self.distanceMatrix.append(arrays)

                lastRow = np.empty(len(self.coordinateArray))
                lastRow.fill(100)

                self.distanceMatrix.append(lastRow)
                self.distanceMatrix = np.array(self.distanceMatrix)
                np.transpose(self.distanceMatrix)
                print len(self.distanceMatrix)
                print len(self.distanceMatrix[0])
#                np.savetxt('coordinateArray.txt', self.coordinateArray) # coordinateArray[0:50,0:50]
#                np.savetxt('np_matrix.txt', self.distanceMatrix) # distanceMatrix[0:50]

        else:
            if rank == 0:
                print ">> less than "+str(size*3)+" solutions, proceeding ..."

                for n in xrange(0,len(self.coordinateArray),1):

                    for m in xrange (n,len(self.coordinateArray),1):
                        # make sure you are not using the same structures against themselves
                        if n == m:
        #                                       # add a "wrong" distance in the matrix to only have half the matrix
                            pass
                        else:

                            # --------------------------------- MODIFY THE FLEXIBLE STRUCTURES FOR THE 1ST ASSEMBLY
                            for structure in self.data.structure_list:
                                if self.params.assembly_style=="flexible" and structure.flexibility != "NA":
                                    len_rigid_dim = 6*(len(self.data.structure_list)-1)
                                    i = 0

                                    deform_coeffs = self.coordinateArray[n][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]

                                    if self.params.mode=="seed":
                                        pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                                        code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                                        target_frame=min_dist.argmin()
                                        coords=structure.flexibility.all_coords[:,target_frame]
                                        coords_reshaped=coords.reshape(len(coords)/3,3)
                                        structure.monomer.set_xyz(coords_reshaped)
                                    else:
                                        coords=structure.monomer.get_xyz()
                                        coords_reshaped=coords.reshape(len(coords)*3)

                                        for n in xrange(0,len(deform_coeffs),1):
                                            coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                                        structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                                    i += len(structure.flexibility.eigenspace_size)
                                else:
                                    structure.monomer.set_xyz(structure.init_coords)

                            # ------------------- CREATING 1ST ASSEMBLY
                            assembly1 = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
                            assembly1.place_all_mobile_structures(self.coordinateArray[n][:len(self.coordinateArray[n])-1])
                            # get the coordinates of all the structures to get the coordinates of Assembly
                            coordinate_of_assembly1_structures = []
                            for structure_index in xrange(0,len(self.data.structure_list), 1):
                                coordinate_of_assembly1_structures.append(assembly1.get_structure_xyz(structure_index))
                            m1 = np.concatenate((coordinate_of_assembly1_structures),axis=0)

                            # --------------------------------- MODIFY THE FLEXIBLE STRUCTURES FOR THE 2ND ASSEMBLY
                            for structure in self.data.structure_list:
                                if self.params.assembly_style=="flexible" and structure.flexibility != "NA":
                                    len_rigid_dim = 6*(len(self.data.structure_list)-1)
                                    i = 0

                                    deform_coeffs = self.coordinateArray[n][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]

                                    if self.params.mode=="seed":
                                        pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                                        code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                                        target_frame=min_dist.argmin()
                                        coords=structure.flexibility.all_coords[:,target_frame]
                                        coords_reshaped=coords.reshape(len(coords)/3,3)
                                        structure.monomer.set_xyz(coords_reshaped)
                                    else:
                                        coords=structure.monomer.get_xyz()
                                        coords_reshaped=coords.reshape(len(coords)*3)

                                        for n in xrange(0,len(deform_coeffs),1):
                                            coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                                        structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                                    i += len(structure.flexibility.eigenspace_size)
                                else:
                                    structure.monomer.set_xyz(structure.init_coords)

                            # ------------------- CREATING 2ND ASSEMBLY
                            assembly2 = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
                            assembly2.place_all_mobile_structures(self.coordinateArray[m][:len(self.coordinateArray[m])-1])
                            # get the coordinates of all the structures to get the coordinates of Assembly
                            coordinate_of_assembly2_structures = []
                            for structure_index in xrange(0,len(self.data.structure_list), 1):
                                coordinate_of_assembly2_structures.append(assembly2.get_structure_xyz(structure_index))
                            m2 = np.concatenate((coordinate_of_assembly2_structures),axis=0)



                            # calculate RMSD between the 2
                            rmsd=self.align(m1,m2) # --> comes from Default.Postprocess.align()
                            self.distanceMatrix[n][m] = rmsd
                            
        if rank == 0:
            np.savetxt('coordinateArray.txt', self.coordinateArray)
            np.savetxt('np_matrix.txt', self.distanceMatrix)
            
            # launch the Clustering and tree drawing module
            app = wx.App(False)
            frame = CnD.MainFrame(None, "Clustering interface",self.OUTPUT_DIRECTORY ,self.params, self.data, self)
            frame.RMSDPanel.computeMatrix()
            if self.params.cluster_threshold == "NA":
                frame.Show()
                app.MainLoop()
            else:
                frame.RMSDPanel.convertCoordsAndExportPDB(self.params.cluster_threshold)
                
                
    def write_pdb(self, centroidArray, average_RMSD_ARRAY):
        iterant = 0 # this iterator is used (in a bad way :( ) to select the right average RMSD value when iterating over the centroids

        clusters_file=open("%s/solutions.dat"%self.OUTPUT_DIRECTORY,"w")

        # writing the tcl file:
        tcl_file = open("%s/assembly.vmd"%self.OUTPUT_DIRECTORY,"w")

        # import the constraint file:
        self.constraint = self.params.constraint.split('.')[0]

        #test if constraint file exists, and function constaint_check is available
        try:
            exec 'import %s as constraint'%(self.constraint)
        except ImportError, e:
            print "ERROR: load of user defined constraint function failed!"
            sys.exit(1)

        try: constraint.constraint_check
        except NameError:
            print 'ERROR: constraint_check function not found'
            
        # HETEROMULTIMER ASSEMBLY
        else:
            print "extracting Complex multimer pdb"
            for n in centroidArray:

                for structure in self.data.structure_list:
                    if self.params.assembly_style=="flexible" and structure.flexibility != "NA":
                        len_rigid_dim = 6*(len(self.data.structure_list)-1) # careful here! the -1 is undecisive
                        i = 0

                        deform_coeffs = self.coordinateArray[n][len_rigid_dim : len_rigid_dim + i + len(structure.flexibility.eigenspace_size) ]

                        if self.params.mode=="seed":
                            pos_eig=structure.flexibility.proj[:,structure.flexibility.centroid]+deform_coeffs
                            code,min_dist=vq(structure.flexibility.proj.transpose(),np.array([pos_eig]))
                            target_frame=min_dist.argmin()
                            coords=structure.flexibility.all_coords[:,target_frame]
                            coords_reshaped=coords.reshape(len(coords)/3,3)
                            structure.monomer.set_xyz(coords_reshaped)
                        else:
                            coords=structure.monomer.get_xyz()
                            coords_reshaped=coords.reshape(len(coords)*3)

                            for n in xrange(0,len(deform_coeffs),1):
                                coords_reshaped+=deform_coeffs[n]*structure.flexibility.eigenvec[:,n]

                            structure.monomer.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                        i += len(structure.flexibility.eigenspace_size)
                    else:
                        structure.monomer.set_xyz(structure.init_coords)
                # ------------------------------ CREATE ASSEMBLY
                print "creating PDB for centroid: "+str(iterant)
                multimer1 = A.AssemblyHeteroMultimer(self.data.structure_list_and_name)
                multimer1.place_all_mobile_structures(self.coordinateArray[n][:len(self.coordinateArray[n])-1])
                # print the pdb file
                multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,iterant))

                # create the constraint:
                measure = constraint.constraint_check(self.data, multimer1)

                # ----------------------------- WRITING SOLUTION.DAT
                l = []
                f = []
                # insert coordinates in the solution.dat file

                f.append("assembly "+str(iterant)+" |")
                for item in self.coordinateArray[n][: len(self.coordinateArray[n])-1]:
                    l.append(item)
                    f.append("%8.3f ")
                    #write constraint values
                f.append("| ")
                for item in measure:
                    l.append(item)
                    f.append("%8.3f ")
                #write fitness
                f.append("| %8.3f")
                l.append(self.coordinateArray[n][-1])
                # write average RMSD OF CLUSTER:
                f.append("| %8.3f\n")
                l.append(average_RMSD_ARRAY[iterant])

                formatting=''.join(f)

                clusters_file.write(formatting%tuple(l))

                # --------------------------- WRITING TCL FILE
                if iterant == 0:
                    tcl_file.write("mol new assembly"+str(iterant)+".pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all \n")

                else:
                    tcl_file.write("mol addfile assembly"+str(iterant)+".pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all \n")

                iterant += 1


        tcl_file.write("mol delrep 0 top \n\
mol representation NewCartoon 0.300000 10.000000 4.100000 0 \n\
mol color Chain \n\
mol selection {all} \n\
mol material Opaque \n\
mol addrep top \n\
mol selupdate 0 top 0 \n\
mol colupdate 0 top 0 \n\
mol scaleminmax top 0 0.000000 0.000000 \n\
mol smoothrep top 0 0 \n\
mol drawframes top 0 {now}")

        clusters_file.close()
        tcl_file.close()
