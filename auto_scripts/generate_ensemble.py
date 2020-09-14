#!/usr/bin/env python

# Author: Lucas Rudden

from Default import Parser as P
from Default import Space as S
from Default import Postprocess as PP

import numpy as np
import sys, os
from sklearn.cluster import KMeans
from copy import deepcopy
import JabberDock as jd
import biobox as bb

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

class Parser(P):
    def __init__(self):
        #protein information - define ligand and receptor
        self.add('receptor','receptor','str',"NA")
        self.add('ligand','ligand','str',"NA")
        # Add in the electron density maps that we need
        self.add('receptor_map','receptor_map','str',"NA")
        self.add('ligand_map','ligand_map','str',"NA")
        # Cutoff for electron density map & distance lists
        self.add('iso_cutoff','iso_cutoff','float',0.5)
        self.add('dist_cutoff','dist_cutoff','float',1.6)
   
        # Clash check for vdw
        self.add('clash_atoms','atoms','array str',["CA", "CB"])
        
        # Number of samples and sampling methods types
        self.add('samples','samples','int',300) 
        self.add('sampling_method','smethod','str', "kmeans")
        self.add('file_name','file_name','str',"file") 

    def check_variables(self):

        #check files existence
        if self.receptor!="NA":
            tmp=os.path.abspath(self.receptor)
            if os.path.isfile(self.receptor)!=1 :
                print("ERROR: receptor building block pdb file %s not found"%self.receptor)
                sys.exit(1)
        if self.ligand!="NA":
            tmp=os.path.abspath(self.ligand)
            if os.path.isfile(self.ligand)!=1 :
                print("ERROR: ligand building block pdb file %s not found"%self.ligand)
                sys.exit(1)
        if self.receptor_map!="NA":
            tmp=os.path.abspath(self.receptor_map)
            if os.path.isfile(self.receptor_map)!=1 :
                print("ERROR: receptor building block dx file %s not found"%self.receptor_map)
                sys.exit(1)
        if self.ligand_map!="NA":
            tmp=os.path.abspath(self.ligand_map)
            if os.path.isfile(self.ligand_map)!=1 :
                print("ERROR: ligand building block dx file %s not found"%self.ligand_map)
                sys.exit(1)



class Data:
    
    def __init__(self,params):
        

        self.params=params

        print("\n> Input DATA:")
        print(">> protein receptor: %s"%params.receptor)
        print(">> ligand: %s"%params.ligand)

        #import protein structures - get receptor and ligand
        print(">> importing PDB files...")
        self.M1 = bb.Molecule()
        self.M2 = bb.Molecule()
        self.M1.import_pdb(params.receptor)
        self.M2.import_pdb(params.ligand)
    
        print(">> importing density files...")
        self.E2 = bb.Density()
        self.E2._import_dx(params.ligand_map) 
        self.E1 = bb.Density()
        self.E1._import_dx(params.receptor_map)

        # Build the jabber maps so we save computation time (essentially just the chosen isosurface from dx)
        self.J1 = jd.Jabber(params.receptor_map, params.iso_cutoff)
        self.J2 = jd.Jabber(params.ligand_map, params.iso_cutoff)                       

        # move the structures to the geometric center (so they're overlapping)
        trans_M1 = deepcopy(self.M1.get_center())
        trans_M2 = deepcopy(self.M2.get_center())

        self.M1.translate(-trans_M1[0], -trans_M1[1], -trans_M1[2])
        self.M2.translate(-trans_M2[0], -trans_M2[1], -trans_M2[2])
        jd.geometry.translate_map(self.E1, -trans_M1[0], -trans_M1[1], -trans_M1[2])
        jd.geometry.translate_map(self.E2, -trans_M2[0], -trans_M2[1], -trans_M2[2])
        self.J1.translate(-trans_M1)
        self.J2.translate(-trans_M2)

        self.E1.write_dx('./models/initial_%s'%(params.receptor_map))
        self.E2.write_dx('./models/initial_%s'%(params.ligand_map))
        self.M1.write_pdb('./models/initial_%s'%(params.receptor))
        self.M2.write_pdb('./models/initial_%s'%(params.ligand))
        self.J1.write_jabber('./models/initial_receptormap.pdb')
        self.J2.write_jabber('./models/initial_ligandmap.pdb')
        ###############################################

class Space(S):
    def __init__(self,params,data):

        # Here we define the search space
        #assign low boundaries
        if params.low_input!="NA" :
            self.low=np.zeros(len(params.low_input))
            for i in xrange(0,len(params.low_input),1):
                self.low[i]=params.low_input[i]

        else:
            print("ERROR: boundaryMin should be defined")
            sys.exit(1) 
        
        #assign high boundaries
        if params.high_input!="NA" :
            self.high=np.zeros(len(params.high_input))
            for i in xrange(0,len(params.high_input),1):
                self.high[i]=params.high_input[i]

        else:
            print("ERROR: boundaryMax should be defined")
            sys.exit(1)
 
 
        ###################################################
 
        #check boundary conditions consistency
        if len(self.low) != len(self.high):
            print 'ERROR: dimensions of min and max boundary conditions are not the same'
            sys.exit(1)
        if (self.low>self.high).any():
            print 'ERROR: a lower boundary condition is greated than a higher one'
            sys.exit(1)
        #define cell size
        self.cell_size=self.high-self.low
                       
        #set boundary type (default is periodic)
        self.boundary_type=np.zeros(len(params.low_input))
        if params.boundary_type!="NA":
            for i in xrange(0,len(params.low_input),1):
                self.boundary_type[i]=params.boundary_type[i]


class Fitness:

    def __init__(self,data,params):
        self.data=data
        self.params=params
        self.m1= self.data.M1.atomselect("*", "*", self.params.atoms)       
        self.COM = self.data.M2.get_center()

    def evaluate(self,num,pos): # num = number of red crosses, pos = specific position in search space (contains parameters in map_manipulation for roto-translations)

        Ptest_pdb = deepcopy(self.data.M2)
        R_pdb = jd.geometry.rotate_pdb(Ptest_pdb, pos[3], pos[4], pos[5], pos[6])
        Ptest_pdb.translate(x = pos[0], y = pos[1], z =pos[2])

        #extract atoms for clash detection, and compute energy
        m2= Ptest_pdb.atomselect("*", "*", self.params.atoms)
        score_check = self.vdw_energy(self.m1, m2)

        if score_check <= 0.0:

            # Import ligand map
            Ptest_J2 = deepcopy(self.data.J2)
            Ptest_J2.rotate(self.COM, R_pdb)
            Ptest_J2.translate(np.array((pos[0], pos[1], pos[2])))
            
            # Assess the score between two maps
            dist, bool_index, min_index = self.data.J1.distance_list(Ptest_J2, cutoff = self.params.dist_cutoff)
            # Negative because we want to minmise the score for POW 
            Sc = - self.data.J1.scoring_function(Ptest_J2, dist, bool_index, min_index)

            if np.isnan(Sc):
                Sc = np.inf
        else:
             Sc = np.inf
        return float(Sc)


    def vdw_energy(self, m1, m2):

        epsilon=1.0
        sigma=2.7
        cutoff=12.0
        energy=0.0
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
        dist=np.array(d)

        # Detect interfacing atoms (atom couples at less than a certain cutoff distance)
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
        for i in xrange(0,len(couples[0]),1):
            d=dist[couples[0,i],couples[1,i]]
            energy+=4.*epsilon*((sigma/d)**9.-(sigma/d)**6.)
            
        return energy


class Postprocess(PP):

    def __init__(self,data,params):
        self.data=data
        self.params=params
        self.fitness=Fitness(data, params)
        self.COM = self.data.M2.get_center()

    def run(self) :

        if rank==0:

            #use superclass method to filter acceptable solutions
            self.log=self.select_solutions(self.params)
            print ">> %s solutions filtered"%len(self.log)

            if len(self.log)==0:
                print ">> no solution found!!!"
                return

            elif len(self.log) > self.params.samples:
                if self.params.smethod == "kmeans":
                    print ">> selecting %s representatives with kmeans..."%self.params.samples
                    kmeans = KMeans(n_clusters=self.params.samples, random_state=0).fit(self.log[:,:-1])
                    points = kmeans.cluster_centers_      

                    # kmeans.labels_ returns the labels for each conformation to point to what cluster (out of samples) it belongs to
                    unique_clusters, unique_index, cluster_count = np.unique(kmeans.labels_, return_index=True, return_counts = True)
                    point_save = self.log[unique_index, :]
                    point_save = np.append(point_save, np.array([cluster_count]).T, axis=1)
                    
                else:
                    print ">> selecting %s random representative..."%self.params.samples
                    pos = np.arange(len(self.log[:, :-1]))
                    selpos = np.random.choice(pos,self.params.samples, False)
                    points = self.log[selpos, :-1]
                    point_save = self.log[selpos]

            else:
                print ">> only %s solutions found..."%len(self.log)
                points = self.log[:,:-1]
                point_save = self.log
            
            point_save[:, 7] *= -1 # convert back so increasing score is better
            # Print the relevent cluster points
            np.savetxt('./models/%s.dat'%(self.params.file_name), point_save)
            
        else:
            points = []
            point_save = []

        
        #broadcast points to everybody    
        comm.Barrier()
        points = comm.bcast(points, root=0)     
        point_save = comm.bcast(point_save, root=0)

        # calculate RMSD all-vs-all
        if rank == 0:
	    #print points
	    print ">> calculating pairwise RMSD of selected samples..."

	    for cnt, pos in enumerate(point_save):
        
                Ptest_pdb = deepcopy(self.data.M2)
                        
                R_pdb = jd.geometry.rotate_pdb(Ptest_pdb, pos[3], pos[4], pos[5], pos[6])
                Ptest_pdb.translate(x = pos[0],y = pos[1],z = pos[2])
                      
                Ptest_pdb.write_pdb("./models/model_%i.pdb"%(cnt))

                # Now we do the density map
                Ptest = deepcopy(self.data.E2)
                R = jd.geometry.rotate_map(Ptest, COM = self.COM, R = R_pdb)
                jd.geometry.translate_map(Ptest, pos[0], pos[1], pos[2])

                Ptest.write_dx("./models/model_%i.dx"%(cnt))
        
        
        comm.Barrier()
