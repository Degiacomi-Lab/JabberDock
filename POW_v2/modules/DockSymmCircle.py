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
from scipy.cluster.vq import *
from copy import deepcopy

from Protein import Protein
import Multimer as M
# import Multimer_CG as MCG
# import flexibility

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

import CCC_simulMap as CCC
import time



class Parser(R):
    def __init__(self):

        #fitness function flags
        self.add('constraint','constraint','str',"NA")
        self.add('target','target','array float',np.array([]))
        self.add('detectClash','detect_clash','str',"on")
        self.add('mixingWeight','mix_weight','float', 0.2)

        #flags for search space and data structure
        self.add('degree','degree','int',-1)
        self.add('radius','radius','float',"NA")

        #postprocessing flags
        self.add('cluster_threshold','cluster_threshold','float',"NA")
        self.add('output_folder','output_folder','str',"result")

        #define style (rigid, flexible)
        self.add('style','style','str',"rigid")

        #flexible assembly flags
        self.add('projection','proj_file','str',"NA")
        self.add('align','align','str',"no")
        self.add('ratio','ratio','float',0.9)
        self.add('trajectory','trajectory','str',"NA")
        self.add('trajSelection','trajselection','array str',np.array(["NA"]))
        self.add('topology','topology','str',"NA")
        self.add('frameDistancePenalty','frame_distance','float',10.0)
        #rigid assembly flags
        self.add('monomer','pdb_file_name','str',"NA")

        #receptor flags
        self.add('receptor','receptor','str',"NA")
        self.add('boundaryMaxReceptor','high_input_rec','array float',"NA")
        self.add('boundaryMinReceptor','low_input_rec','array float',"NA")
        self.add('z_padding','pad','int',5.0)
        
        #constraint flags and boudary 
        self.add('constrained_opt','constrained_opt','str','off')
        self.add('targetLow','targetLow','array float',np.array([]))
        self.add('targetHigh','targetHigh','array float',np.array([]))
        
        
        #forcefield flag
        self.add('AllAtom','AllAtom','str','off')
        self.add('topFile','topFile','str','NA')
        self.add('prmFile','prmFile','str','NA')
        self.add('CoarseGrain', 'CoarseGrain', 'str', 'off')
        self.add('cg_file', 'cg_file', 'str', 'NA')
        
        # density map docking flag
        self.add('density_map_docking','map_dock_OnOff', 'str', 'no')
        self.add('density_map','density_map', 'str', 'NA')
        self.add('map_resolution', 'map_resolution', 'float', 10.0 )
        self.add('map_docking_weight', 'map_docking_weight', 'float', 20.0 )
        
        
        # flag for direct rmsd comparison: 07/05/14
        self.add('oriPDB','oriPDB', 'str', 'NA')
        self.add('directRMSDcomp','directRMSDcomp', 'str', 'off')
        
        

    def check_variables(self):

        #needed, set flag of clash detection
        if self.detect_clash == 'on' :
            self.detect_clash=True
        elif self.detect_clash== 'off':
            self.detect_clash=False
        else:
            print 'ERROR: detectClash should either be "on" or "off"!'
            sys.exit(1)

        #check if fitness function's energy vs geometry mixing weight is sane
        if self.mix_weight>1 or self.mix_weight<0:
            print 'ERROR: mixingWeight should be a float within 0 and 1!'
            sys.exit(1)

        #check if number of monomers is provided
        if self.degree==-1 :
            print 'ERROR: number of monomer not provided!'
            sys.exit(1)

        #check if number of target measures are provided
        if len(self.target) == 0 and not self.constrained_opt:
            print 'ERROR: target measures not specified!'
            sys.exit(1)
            
        #check constraint flag
        if self.constrained_opt == 'on' :
            self.constrained_opt=True
        elif self.constrained_opt== 'off':
            self.constrained_opt=False
        else:
            print 'ERROR: constraint should either be "on" or "off"!'
            sys.exit(1)
            
            
        #check if number of target measures are provided
        if (len(self.targetLow) == 0 or len(self.targetHigh) == 0) and self.constrained_opt:
            print self.constrained_opt
            print 'ERROR: target low and high measures not specified with constrained optimization!'
            sys.exit(1)

        #verify existence of constraint file (in search path)
        if self.constraint=="NA":
            print "ERROR: constraint file not provided"
            sys.exit(1)
        if os.path.isfile(self.constraint)!=1:
            print "ERROR: constraint file %s not found"%self.constraint
            sys.exit(1)

        #add path of constraint file to search space
        a,b=os.path.split(self.constraint)
        self.constraint=b
        sys.path.append(a)

        #test if constraint file exists, and function constaint_check is available
        try:
            exec 'import %s as constraint'%(self.constraint.split('.')[0])
        except:
            print "ERROR: load of user defined constraint function failed!"
            print "ERROR: awww, there may be a syntax error in there..."
            sys.exit(1)
        try:
            constraint.constraint_check
        except AttributeError:
            print 'ERROR: constraint_check function not found in file %s'%self.constraint
            sys.exit(1)


        #check assembly style flag
        if self.style!="rigid" and self.style!="flexible":
            print 'ERROR: style should either be "rigid" or "flexible"!'
            sys.exit(1)
            
        ###RIGID ASSEMBLY FLAGS###
        if self.pdb_file_name=="NA" and self.style=="rigid":
            print "ERROR: in rigid style (default), a monomer should be provided!"
            sys.exit(1)
        elif self.pdb_file_name!="NA" and self.style!="rigid":
            print "ERROR: %s style requested"%self.style
            print "ERROR: usage of monomer flag requires rigid style!"
            sys.exit(1)
        elif os.path.isfile(self.pdb_file_name)!=1 and self.style=="rigid":
            print "ERROR: pdb file %s not found"%self.pdb_file_name
            sys.exit(1)


        ###FLEXIBLE ASSEMBLY FLAGS###
        #verify existence of topology
        if self.topology!="NA":
            if self.style=="flexible":
                tmp=os.path.abspath(self.topology)
                self.topology=tmp
            else:
                print "ERROR: %s style requested"%self.style
                print "ERROR: usage of topology flag requires flexible style!"
                sys.exit(1)
        if os.path.isfile(self.topology)!=1 and self.style=="flexible":
            print "ERROR: topology %s not found"%self.topology
            sys.exit(1)

        #verify existence of traj
        if self.trajectory!="NA":
            if self.style=="flexible":
                tmp=os.path.abspath(self.trajectory)
                self.trajectory=tmp
            else:
                print "ERROR: %s style requested"%self.style
                print "ERROR: usage of trajectory flag requires flexible style!"
                sys.exit(1)
        if os.path.isfile(self.trajectory)!=1 and self.style=="flexible":
            print "ERROR: trajectory %s not found"%self.trajectory
            sys.exit(1)

        #verify existence of projections file
        if self.proj_file!="NA":
            if self.style=="flexible":
                tmp=os.path.abspath(self.proj_file)
                self.proj_file=tmp
            else:
                print "ERROR: %s style requested"%self.style
                print "ERROR: usage of projection flag requires flexible style!"
                sys.exit(1)
        if self.proj_file!="NA" and os.path.isfile(self.proj_file)!=1 and self.style=="flexible":
            print "ERROR: projections file %s not found"%self.proj_file
            sys.exit(1)

        #verify format trajectory selection (if needed)
        if self.trajselection[0]!="NA" and self.style!="flexible":
            print "ERROR: %s style requested"%self.style
            print "ERROR: usage of trajSelection flag requires flexible style!"
            sys.exit(1)
        if self.style=="flexible" and self.trajselection[0]=="NA":
            tmp="protein"
        else:
            tmp =' '.join(self.trajselection[0:len(self.trajselection)])
        self.trajselection=tmp


        #RECEPTOR FLAGS
        if self.receptor!="NA":

            #verify existence of receptor
            tmp=os.path.abspath(self.receptor)
            self.receptor=tmp
            if os.path.isfile(self.receptor)!=1:
                print "ERROR: receptor file %s not found"%self.receptor
                sys.exit(1)

            #check receptor boundaries consistency
            if len(self.high_input_rec)!=len(self.low_input_rec):
                print "ERROR: boundaryMinReceptor and boundaryMaxReceptor should have the same length!"
                sys.exit(1)

            #if len(self.high_input_rec)!=len(self.boundary_type_rec) and self.boundary_type_rec!="NA":
            #    print "ERROR: boundaryTypeReceptor length inconsistent with low and high boundary length!"
            #    print 'ERROR: %s dimensions given, but %s needed!'%(len(self.boundary_type),len(self.low_input))
            #    sys.exit(1)

            if self.high_input_rec!="NA" and (self.low_input_rec>self.high_input_rec).any():
                print 'ERROR: a receptor lower boundary condition is greated than a higher one!'
                print self.low_input_rec
                print self.high_input_rec
                sys.exit(1)
                
        #FORCE-FIELD FLAG
        if self.AllAtom == 'off':
            self.AllAtom = False
        else :
            self.AllAtom = True
            # check if the prm and top files exist
            if os.path.isfile(self.prmFile)!=1:
                print "ERROR: paramater file %s not found"%self.prmFile
                sys.exit(1)
            if os.path.isfile(self.topFile)!=1:
                print "ERROR: topology file %s not found"%self.self.topFile
                sys.exit(1)
                
        if self.CoarseGrain == 'off':
            self.CoarseGrain = False
        
        else :
            self.CoarseGrain = True
            # check if the prm and top files exist
            if os.path.isfile(self.cg_file)!=1:
                print "ERROR: Martini converter python file %s not found"%self.cg_file
                sys.exit(1)
                
        # check whether density flag is on or off
        if self.map_dock_OnOff == "on" :
            if self.density_map != "NA":
                self.map_dock_OnOff = True
            else:
                print 'ERROR: if density map flag is on, a density map should be present in the directory'
                sys.exit(1)
            # check on the parameters:
            if self.map_resolution < 0:
                print 'ERROR: negative resolutions do not exist ...'
                sys.exit(1)
            if self.map_docking_weight < 0:
                print 'ERROR: negative weights not permitted here ...'
                sys.exit(1)
            
            print ">> Electron Density Map Docking Mode"
        elif self.map_dock_OnOff == "off":
            self.map_dock_OnOff = False
        else:
            print 'ERROR: density_map_docking should be either "on" or "off"'
            sys.exit(1)
            
        
            
        
            
        # CHECKING whether to do direct rmsd calculations:
        # check whether density flag is on or off
        #FORCE-FIELD FLAG
        if self.directRMSDcomp == 'off':
            self.directRMSDcomp = False
        else :
            self.directRMSDcomp = True
            # check if the prm and top files exist
            if os.path.isfile(self.oriPDB)!=1:
                print "ERROR: paramater file %s not found"%self.prmFile
                sys.exit(1)
            


class Data:

    def __init__(self,params):

        self.index=[]
        self.index_receptor=[]

        self.eigenspace_size=""

        if params.style=="flexible":
            print ">> flexible docking requested, launching PCA..."
            #check align flag
            if params.align!="yes" and params.align!="no":
                print "ERROR: align flag should be set either to \"yes\" or \"no\"!"
                sys.exit(1)

            #load trajectory
            try:
                self.traj=flexibility.loader(params.trajectory, params.topology, params.align, params.trajselection)
            except ImportError, e:
                print "ERROR: Trajectory loading failed, aborting!"
                sys.exit(1)

            #if projection on eigenvectors has been already performed, use the user provided file
            if params.proj_file!="NA":
                if os.path.isfile(params.proj_file)!=1 :
                    print "ERROR: projections file not found!"
                    sys.exit(1)
                else:
                    print ">> loading frames projections file..."
                    p_file = open(params.proj_file, 'r+')
                    p_list=[]
                    line = p_file.readline()
                    while line:
                        w = line.split()
                        p_list.append(w)
                        line = p_file.readline()
                    self.proj=np.array(p_list).transpose().astype(float)
                    p_file.close()
                    if self.traj.shape[1]!=self.proj.shape[1]:
                        print "ERROR: trajectory has %s frames, projection has %s, aborting!"%(self.traj.shape[1],self.proj.shape[1])
                        sys.exit(1)
            #else, peform PCA on trajectory
            else:
                try:
                    self.proj,e=flexibility.pca(params.trajectory,params.topology, params.ratio)
                except ImportError, e:
                    print "ERROR: PCA failed, aborting!"
                    sys.exit(1)


            #compute projections centroid
            try:
                self.centroid=flexibility.clusterize(self.proj)

            except ImportError, e:
                print "ERROR: Clustering failed, aborting!"
                sys.exit(1)

            #suppose that data assigned to cluster can be fitted with a gaussian.
            #interval will be set at 95% of obtained distribution
            p=[]
            for i in xrange(0,len(self.proj),1):
                p.append(np.std(self.proj[i])*2)

            self.eigenspace_size=np.array(p)
            print ">> Flexibility interval in 95% confidence..."
            print "   "+str(self.eigenspace_size)

            #PCA automatically generates a reference pdb file called "protein.pdb"
            #check its existence and load its name in local params object!
            params.pdb_file_name="protein.pdb"
            if os.path.isfile(params.pdb_file_name)!=1 :
                print "ERROR: pdb file %s not found"%params.pdb_file_name
                sys.exit(1)

        #load monomeric structure (both for rigid )
        self.structure = Protein()
        try:
            self.structure.import_pdb(params.pdb_file_name)
        except IOError:
            sys.exit(1)

        #load receptor structure, if provided, and center it to the origin
        self.receptor_structure=[]
        if params.receptor!="NA":
            print ">> loading receptor, and centering it to origin..."
            self.receptor_structure = Protein()
            try:
                self.receptor_structure.import_pdb(params.receptor)
            except IOError,e:
                print e
                print "Man, a loading error!"
                sys.exit(1)

            coords=self.receptor_structure.get_xyz()
            self.receptor_structure.set_xyz(coords-np.mean(coords,axis=0))


        #if clash avoidance is requested, compute indexes of relevant atom
        #this call currently just selects the CA atoms
        if params.detect_clash and not params.CoarseGrain:
            #self.index=self.get_index(["CA","CB"])
            #self.index=self.get_index(["CA"])
            self.get_index(["CA"])
            
        # in case we re doing ALL ATOM docking:
        if params.AllAtom:
            self.data_list = np.array(self.structure.mapping(self.structure.data))[:,1:5]
            self.AtomVDWprm, self.ResidueAtomCharge = self.extract_AllAtom_info(params)
            
        # in case of COARSE GRAINED modeling:
        if params.CoarseGrain:
            import subprocess
            
            # create the CG model of the multimer:
            #subprocess.call('python '+params.cg_file+' -f '+params.pdb_file_name+' -x cg_model.pdb -p backbone -ff martini22', shell=True)
            self.cg_structure = Protein()
#             self.cg_structure.import_pdb('cg_model.pdb')
            self.cg_structure.import_pdb('CG_protein_0.pdb')
            print ">> done making cg model of monomer ..."
            
        # if the density map docking is on load the structure into data:
        if params.map_dock_OnOff:
            self.density_map_fileName = params.density_map
            
        # if the direct RMSD comparison flag is on setting up data
        if params.directRMSDcomp:
            self.rmsdAtoms = self.extract_data_for_directRMSD(params)
            self.Fevals     = 0


    def get_index(self,atoms=["CA","CB"]):

        #generate a dummy multimer
        multimer = M.Multimer(self.structure)
        multimer.create_multimer(3,10,np.array([0.0,0.0,0.0]))

        #extract the indexes where atoms of interest are located
        for aname in atoms:
            [m,p_index]=multimer.atomselect(1,"*","*",aname,True)
            for i in p_index:
                self.index.append(i)

#             if self.receptor_structure!=[]:
#                 [m,p_index]=self.receptor_structure.atomselect("*","*",aname,True)
#                 for i in p_index:
#                     self.index_receptor.append(i)
        #return self.index
        
    def extract_AllAtom_info(self,params):
        '''This function parses the .prm and .top files of the corresponding force field and extract the charges, Emin and Rmin of each atom 
        residue to compute the NON BONDED interactions between proteins monomers'''
        
        
        # ----------------------------------------- PARSING PRM FILE CHARMM -----------------------------------------

        prm_file = open(params.prmFile)
        
        # first parse the topology file to extract atom Epsilon and rmin
        writing_signal = 0 # will turn to 1 as soon as we encounter the keyword "NONBONDED" in the prm file!
        prmAtom = {} # -> {atomName: [Emin, rmin], ...}
        
        
        while True:
            line = prm_file.readline()
            if line:
                # activate and deactivate the writing signal depending on the header encountered
                splittedLine = line.split(" ")
                if splittedLine[0] == "NONBONDED": # approx at line 3436
                    writing_signal = 1
                elif splittedLine[0] == "HBOND": # terminate when get to this flag
                    writing_signal = 0
                    
                #when you encounter the description of the non bonded interaction, start the recording
                elif writing_signal and splittedLine[0] != "cutnb" and line[0] != "!" \
                and splittedLine[0] != "\t\t!" and splittedLine[0] != "" and splittedLine[0] != "\n":
                    counter=0 # to make sure you only take the Emin and rmin!
                    # record the information stored in each line!
                    prmAtom[splittedLine[0]] = []
                    for i in xrange(0,len(splittedLine),1):
                        if splittedLine[i] != "" and i != 0 and counter < 3 :
                            if counter == 0:
                                counter +=1
                            else:
                                prmAtom[splittedLine[0]].append(float(splittedLine[i]))
                                counter += 1
            else:
                break
        
        
        # ----------------------------------------- PARSING TOP FILE CHARMM -----------------------------------------
                 
        top_file = open(params.topFile)
        
        writing_signal = 0 # same usage as above
        topResidue = {} # structure -> {"topResidue1": { atomName1(PDB): [atomName1(.prm), charge!] , atomName2(PDB): [atomName2(.prm), charge!], ... }, "topResidue2":{}, ... } 
        
        while True:
            line = top_file.readline()
            if line:
                splittedLine = line.split(" ")
                
                if splittedLine[0] == "RESI":
                    writing_signal = 1
                    currentResidue = splittedLine[1]
                    topResidue[currentResidue] = {}
                elif splittedLine[0] == "PRES" and splittedLine[1] == "GLYP":
                    break
                
                if writing_signal and splittedLine[0] == "ATOM":
                    topResidue[currentResidue][splittedLine[1]] = []
                    
                    # add the prm nomenclature as well as the charge to the array:
                    counter=0
                    for word in splittedLine:
                        if word != "" and word != "ATOM" and counter < 3:
                            if counter == 0:
                                counter +=1
                            else:
                                wordstripped=word.replace("\n","")
                                topResidue[currentResidue][splittedLine[1]].append(wordstripped)
                                counter += 1
                    
                    # render float the charge value which is string      
                    topResidue[currentResidue][splittedLine[1]][-1] = float(topResidue[currentResidue][splittedLine[1]][-1])
                    counter += 1
            else:
                break
        
        return prmAtom, topResidue
        
    def extract_data_for_directRMSD(self,params):
        
        # ------------- Ref Protein
        ref = Protein()
        ref.import_pdb(params.oriPDB)
        data_list = ref.mapping(ref.data)
        # create dictionary composed of :
        # _> {0: 'N_MET_A_3', 1: 'CA_MET_A_3', ...}
        self.ref_pdb_dictio = {}
        for i in xrange(0,len(data_list),1):
            inf = data_list[i]
            self.ref_pdb_dictio[i] = inf[1]+"_"+inf[2]+"_"+inf[3]+\
                                       "_"+str(inf[4]) 
                                       
        # ------------- Model Protein
        # then extract index of atoms to be aligned for both ref and model
        # first create dummy multimer
        
        multimer      = M.Multimer(self.structure)
        multimer.create_multimer(params.degree,10,np.array([0.0,0.0,0.0]))
        
        # --------- tests
#         exampledata = [67.4128881891, 277.013912548, 256.203841598, -0.767253828527]
#         multimer.create_multimer(params.degree,exampledata[3],np.array([exampledata[0],exampledata[1],exampledata[2]]))
#         multimer.write_PDB('assembly1.pdb')
        # ---------
        
        multimer_list = []
        
        chain_converter = ('A','B','C','D','E','F','G','H','I','J','K','L','M',\
                           'N','O','P','Q','R','S','T','U','V','W','X','Y','Z')
        
        
        xyz_mult_list = []
        for j in xrange(0,params.degree,1):
            multimer.pdb.set_xyz(multimer.multimer[j])

            #map intergers to characters from input data (default: all the protein)
            tmp_list=multimer.pdb.mapping(multimer.pdb.data)
            temp = np.array(tmp_list)
            temp[:,3] = chain_converter[j]
            multimer_list += temp
            xyz = [[x[5],x[6],x[7]] for x in tmp_list]
            xyz_mult_list += xyz 
        
        self.model_pdb_dictio = {}
        
        for i in xrange(0,len(multimer_list),1):
            inf = multimer_list[i]
            self.model_pdb_dictio[inf[1]+"_"+inf[2]+"_"+inf[3]+\
                                       "_"+str(inf[4])] = i
        
        # ------------- Atom extraction
        # check which atoms the ref and model have in common
        self.RMSDcommonAtoms = {}
        self.RMSDcommonAtoms['ref'] = []
        self.RMSDcommonAtoms['mod'] = []
        
        for i in xrange(0,len(self.ref_pdb_dictio),1):
            atomName = self.ref_pdb_dictio[i].split('_')[0]
            if atomName == 'CA':
                if self.ref_pdb_dictio[i] in self.model_pdb_dictio.keys():
                    self.RMSDcommonAtoms['ref'].append(i)
                    self.RMSDcommonAtoms['mod'].append(\
                        self.model_pdb_dictio[self.ref_pdb_dictio[i]] )
        
        # example
        xyz_data_list = [[x[5],x[6],x[7]] for x in data_list]
        npData_list     = np.array(xyz_data_list)
        
        npMult_list = np.array(xyz_mult_list)
        
        
        self.xyz_ref = npData_list[self.RMSDcommonAtoms['ref']]
        xyz_mod = npMult_list[self.RMSDcommonAtoms['mod']]
        
#         print self.get_RMSD(self.xyz_ref, xyz_mod)
           
#         ccc
                
            
    def get_RMSD(self,m1, m2):

        L = len(m1)

        #centering protein
        COM1 = np.sum(m1,axis=0) / float(L)
        COM2 = np.sum(m2,axis=0) / float(L)
        m1 -= COM1
        m2 -= COM2

        E0 = np.sum( np.sum(m1*m1,axis=0),axis=0) + np.sum( np.sum(m2*m2,axis=0),axis=0)

        #This beautiful step provides the answer. V and Wt are the orthonormal
        # bases that when multiplied by each other give us the rotation matrix, U.
        # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
        V, S, Wt = np.linalg.svd( np.dot( np.transpose(m2), m1))

        reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:,-1] = -V[:,-1]

        RMSD = E0 - (2.0 * sum(S))
        RMSD = np.sqrt(abs(RMSD / L))

        return RMSD
        
        

class Space(S):
    # n dim. space: 3 monomer rotations, 1 monomer radial translation, 1 multimer z translation, 3 receptor rotations, n monomer eigenvectors
    def __init__(self,params,data):

        #add receptor dimensions if needed
        rec_dim=0
        if params.receptor!="NA":
            rec_dim=4

        self.low=np.zeros(4+rec_dim+len(data.eigenspace_size))
        self.high=np.zeros(4+rec_dim+len(data.eigenspace_size))
        self.cell_size=np.zeros(4+rec_dim+len(data.eigenspace_size))
        self.boundary_type=np.zeros(4+rec_dim+len(data.eigenspace_size))

        if len(params.high_input)!=len(params.low_input):
            print "ERROR: seems that low and high boundaries have different length..."
            sys.exit(1)

        #assign monomer low boundaries
        if params.low_input!="NA" :
            if len(params.low_input)==3 or len(params.low_input)==4:
                for i in xrange(0,len(params.low_input),1):
                    self.low[i]=params.low_input[i]
            else:
                print "ERROR: boundaryMin should contain within 3 or 4 values"
                sys.exit(1)
        else:
            print "WARNING: boundaryMin undefined, using default values for monomer rotation"

        #assign monomer high boundaries
        if params.high_input!="NA" :
            if len(params.high_input)==3 or len(params.high_input)==4:
                for i in xrange(0,len(params.high_input),1):
                    self.high[i]=params.high_input[i]
            else:
                print "ERROR: boundaryMax should contain 3 or 4 values"
                sys.exit(1)
        else:
            print "WARNING: boundaryMin undefined, using default values for monomer rotation"
            self.high[0]=360.0
            self.high[1]=180.0
            self.high[2]=360.0

        #if needed assign radius boundaries
        if params.radius!="NA":
            if params.high_input=="NA":
                self.low[3]=params.radius
                self.high[3]=params.radius
            else:
                print "ERROR: radius keyword used when also using boundary condition list!"
                sys.exit(1)
        if params.radius=="NA" and params.high_input=="NA":
            print "ERROR: translation undefined boundary or radius keywords needed!"
            sys.exit(1)


        #if needed, add receptor dimensions
        if params.receptor!="NA":

                #multimer z-axis translation
            if params.high_input_rec=="NA":
                d=data.structure.get_xyz()
                self.low[4]=d[2].min()+params.pad
                self.high[4]=d[2].max()+params.pad
            else:
                self.low[4]=params.low_input_rec[0]
                self.high[4]=params.high_input_rec[0]

                #receptor x,y,z axis rotations
                if len(params.high_input_rec)!=4:
                    print "WARNING: Receptor rotations undefined in boundaryMinReceptor, using default values"
                    self.low[5]=0.0
                    self.low[6]=0.0
                    self.low[7]=-360/(2*params.degree)
                    self.high[5]=0
                    self.high[6]=0
                    self.high[7]=360/(2*params.degree)
                else:
                    self.low[5]=params.low_input_rec[1]
                    self.high[5]=params.high_input_rec[1]
                    self.low[6]=params.low_input_rec[2]
                    self.high[6]=params.high_input_rec[2]
                    self.low[7]=params.low_input_rec[3]
                    self.high[7]=params.high_input_rec[3]


        #add eigenvector fluctuations in search space
        #for i in xrange(0,len(data.eigenspace_size),1):
        #    self.low[4+rec_dim+i]=-data.eigenspace_size[i]
        #    self.high[4+rec_dim+i]=data.eigenspace_size[i]

        #print "PROJ SHAPE: %s %s"%(data.proj.shape[0],data.proj.shape[1])
        #for i in xrange(0,len(data.eigenspace_size),1):
        #    self.low[4+rec_dim+i]=-data.eigenspace_size[i]
        #    self.high[4+rec_dim+i]=data.eigenspace_size[i]

        #print "PROJ SHAPE: %s %s"%(data.proj.shape[0],data.proj.shape[1])
        for i in xrange(0,len(data.eigenspace_size),1):
            self.low[4+rec_dim+i]=data.proj[i,:].min()
            self.high[4+rec_dim+i]=data.proj[i,:].max()
        #final check for all boundary conditions consistency (cause we're paranoid)
        if (self.low>self.high).any():
            print 'ERROR: a lower boundary condition is greater than a higher one'
            sys.exit(1)

        #define cell size
        self.cell_size=self.high-self.low

        #set boundary type
        if params.boundary_type[0]=="NA":
            for i in xrange(0,len(self.low),1):
                self.boundary_type[i]=0
        elif params.boundary_type[0]!="NA" and len(params.boundary_type)!=len(self.low):
            print 'ERROR: boundaries type inconsistent with system dimensions'
            print 'ERROR: %s dimensions given, but %s needed!'%(len(params.boundary_type),len(self.low))
            sys.exit(1)
        else:
            for i in xrange(0,len(self.low),1):
                self.boundary_type[i]=params.boundary_type[i]


class Fitness:

    def __init__(self,data,params):

        self.params=params
        self.style=params.style
        self.clash=params.detect_clash
        if self.params.constrained_opt:
            self.targetLow=params.targetLow
            self.targetHigh=params.targetHigh
            self.target=params.target 
        else:
            self.target=params.target
        self.map_docking_flag = params.map_dock_OnOff # do so because you want to pass this var to the function evaluate below

        # loading the reference/experimental density map file if flag is on
        if self.map_docking_flag:
            self.density_map_fileName = params.density_map
            

        self.constraint=params.constraint.split('.')[0]
        #test if constraint file exists, and function constaint_check is available
        try:
            exec 'import %s as constraint'%(self.constraint)
        except ImportError, e:
            print "ERROR: load of user defined constraint function failed!"
            sys.exit(1)

        try:
            constraint.constraint_check
        except AttributeError:
            print 'ERROR: constraint_check function not found'
            sys.exit(1)

        #data to manipulate
        self.data=data
        self.degree=params.degree

        #receptor space size
        self.rec_dim=0
        if params.receptor!="NA":
            self.rec_dim=4

        #energy vs geometry mixing weight
        self.c1=params.mix_weight
        
        # parameters for density map docking:
        self.c2 = self.params.map_docking_weight # try make vary from 1 to 100
        self.resol = self.params.map_resolution
        
        # loading the reference/experimental density map file if flag is on
        


    def evaluate(self,num,pos):

        #import user provided constraint file and class needed to construct multimer
        exec 'import %s as constraint'%(self.constraint)
        import Multimer as M

        frame_penalty=0.0

        if self.style=="flexible":
        #pick monomeric structure from database
            deform_coeffs=pos[(4+self.rec_dim):len(pos)]
            #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
            code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
            #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
            target_frame=min_dist.argmin()
            coords=self.data.traj[:,target_frame]
            coords_reshaped=coords.reshape(len(coords)/3,3)
            self.data.structure.set_xyz(coords_reshaped)

        #penalize solutions too far from the sampled space
        #if min_dist.min()>self.params.frame_distance:
        #    frame_penalty=min_dist.min()-self.params.frame_distance
        #print "penalty %s : %s"%(num,frame_penalty)

        #assemble multimer
        if self.params.CoarseGrain:
            self.multimer = MCG.Multimer(self.data.structure, self.data.cg_structure)
            self.multimer.create_multimer(self.degree,pos[3],pos[0:3])
#             self.multimers.write_PDB("exampleAA.pdb")
#             self.multimers.write_CG("exampleCG.pdb")
#             import subprocess
#             subprocess.call('python martinize.py -f exampleAA.pdb -x cg_from_AA.pdb -p backbone -ff martini22', shell=True)
            
        else:
            if self.params.directRMSDcomp:
                self.data.structure = Protein()
                self.data.structure.import_pdb(self.params.pdb_file_name)
            self.multimer = M.Multimer(self.data.structure)
            self.multimer.create_multimer(self.degree,pos[3],pos[0:3])
        
        #------- this part is to write the pdb of interest from coordinates -------
#         singlePos = []
#         for i in pos:
#             singlePos.append(str(i)) 
#         concatPos = "_".join(singlePos)
#         print ">> writing pdb file from coords ..."
# #         self.multimer.write_PDB(concatPos+".pdb")
#         self.multimer.write_PDB("assembly_coordReader.pdb")
        

        # ------ Calculating RMSD on the fly
        if self.params.directRMSDcomp:
            # extract correct atoms for rmsd comparison
            tmp     = self.multimer.get_multimer_uxyz()
            mod_all = np.concatenate((tmp), axis=0)
            xyz_mod = mod_all[self.data.RMSDcommonAtoms['mod']]
                        
            rmsd  =  self.data.get_RMSD(self.data.xyz_ref, xyz_mod)
            
            # write the rmsds in the appropriate file:
            if self.data.Fevals == 0:
                self.rmsd_file = open("rmsd.txt","w")
            
            if self.data.Fevals < self.params.maxFevals:
                self.rmsd_file.write(str(rmsd)+"\n")
                
            self.data.Fevals += 1
            
            if self.data.Fevals > self.params.maxFevals:
                self.rmsd_file.close()
            
#             print rmsd
#             '> computed rmsd and structure file, exiting ...'
#             sys.exit(1)
            
            
            
            
        

        measure=[]
        if self.data.receptor_structure!=[]:
            #shift multimer on z axis
            self.multimer.z_to_origin()
            self.multimer.z_shift(pos[4])

            #rotate receptor (backing up original position, and putting it back when measure is done)
            coords=self.data.receptor_structure.get_xyz()
            self.data.receptor_structure.rotation(pos[5],pos[6],pos[7])
            measure = constraint.constraint_check(self.multimer,self.data.receptor_structure)
            self.data.receptor_structure.set_xyz(coords)
        else:
            try:
                measure = constraint.constraint_check(self.multimer)
            except:
                print "ERROR: syntax error in constraint file!"
                raise
        
        if len(measure) != len(self.target) :
            print 'ERROR: constraint file produced %s measures, but %s target measures are provided!'%(len(measure),len(self.target))
            sys.exit(1)
        
        #measure distance within target values and obtained values
        diff=self.target-np.array(measure)
        distance=np.sqrt(np.dot(diff,diff))
        
        measureArray = list(measure)
        
        # --- Constraint optimization
        if self.params.constrained_opt:
            
            #detect clashes
            if self.clash:
                # in the following cases we want to compute the All atom non bonded only when boundaries are converged, when num is 1!
                if self.params.AllAtom and num == 1:
                    energy = self.AllAtomNonBonded()
                elif self.params.AllAtom and num == 0:
                    energy = np.inf
                elif self.params.CoarseGrain :           
                    energy = self.CGNonBonded(pos)
#                     energy = self.interface_vdw()
                elif self.map_docking_flag:
                    # create the pbd file to be transformed into the density map
                    if self.params.updateSimultaneous:
                        self.multimer.create_PDB_for_density_map(rank, self.params.pdb_file_name)
                        energy = 1 - CCC.get_DMAP_fit( str(rank), self.resol ) 
                        print "> 1 - ccc: "+str(energy)
                
                    else:
                        if num == 1:
                            self.multimer.create_PDB_for_density_map(rank, self.params.pdb_file_name)
                            energy = 1 - CCC.get_DMAP_fit( str(rank), self.resol ) 
                            print "> 1 - ccc: "+str(energy)
                        else:
                            energy=self.interface_vdw()
                else:
                    energy=self.interface_vdw()
                
                # this is the list of the all measure to be sent to CMAVIE and DE 
                measureArray.append(energy)
                
                # this is the single single fitness function used before
                fitness = self.c1*energy+(1-self.c1)*distance+frame_penalty
                
                return measureArray, fitness  
            
           
            else:
                return distance+frame_penalty
          
        # --- Single fitness function optimization  
        else:
            
            #detect clashes
            if self.clash:
                if self.params.AllAtom:
#                     energy = self.AllAtomNonBonded()
                    energy = self.AllAtomNonVW()
                    
                elif self.params.CoarseGrain: 
                    energy = self.CGNonBonded(pos)
                
                elif self.map_docking_flag:
                    # create the pbd file to be transformed into the density map
                    start = time.time() # TIME
                    self.multimer.create_PDB_for_density_map(rank, self.params.pdb_file_name)
                    energy = 1 - CCC.get_DMAP_fit( str(rank), self.resol ) 
                    print "> 1 - ccc: "+str(energy)
                    laps = time.time()-start # TIME
                    print "> process took per assembly: "+str(laps)
                    print "> on 20k evals (hrs): "+str(((20000.0*laps)/60.0)/60)
                            
                else:
                    energy=self.interface_vdw() 
                    
                fitness_score = self.c1*energy+(1-self.c1)*distance+frame_penalty
                
                # this was added for benchmark purposes
#                 measureArray.append(energy)
#                 if self.params.directRMSDcomp:
#                     return  fitness_score, measureArray, rmsd
                    
                if self.params.directRMSDcomp:
                    return  fitness_score
                else: 
                    return  fitness_score
                
            else:
                return distance+frame_penalty
    

    def interface_vdw(self): 

        epsilon=1.0
        sigma=4.7
        cutoff=12.0

        energy=0

        #CLASH WITH NEIGHBORING MONOMER
        #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms

        m1=self.multimer.get_multimer_uxyz()[0][self.data.index]
        m2=self.multimer.get_multimer_uxyz()[1][self.data.index]

        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

        dist=np.array(d)

        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues

        for i in xrange(0,len(couples[0]),1):
            d=dist[couples[0,i],couples[1,i]]
            energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6) ## was 9

        #CLASH WITH RECEPTOR (if needed)
        if self.data.receptor_structure!=[]:
            m1=self.data.receptor_structure.get_xyz()[self.data.index_receptor]
            #extract distances of every atom from all the others
            d=[]
            for i in xrange(0,len(m1),1):
                d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
            dist=np.array(d)
            #detect interfacing atoms (atom couples at less than a certain cutoff distance
            couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
            

            for i in xrange(0,len(couples[0]),1):
                d=dist[couples[0,i],couples[1,i]]
                energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
        
        return energy

    def clashes(self):
        epsilon=1.0
        sigma=4.7
        cutoff=2.5

        energy=0

        #CLASH WITH NEIGHBORING MONOMER
        #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms

        m1=self.multimer.get_multimer_uxyz()[0]
        m2=self.multimer.get_multimer_uxyz()[1]

        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

        dist=np.array(d)

        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
        
        # create an array to store the clashes and their respective energy
        infoData = [] # -> [[atom1, atom2, clashenergy],..]
        
        datalist = np.array(self.data.structure.mapping(self.data.structure.data))[:,1:5]
        
        for i in xrange(0,len(couples[0]),1):
            d=dist[couples[0,i],couples[1,i]]
            energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
            nrj = 4*epsilon*((sigma/d)**9-(sigma/d)**6)
            infoAtom1 = datalist[couples[0,i]]
            infoAtom2 = datalist[couples[1,i]]
            infoData.append([couples[0,i],infoAtom1[0],infoAtom1[1],infoAtom1[2],infoAtom1[3],couples[1,i],infoAtom2[0],infoAtom2[1],infoAtom2[2],infoAtom2[3],nrj])
        
        
        return energy, infoData
        
    
    def AllAtomNonBonded(self):
        cutoff = 12.0 # the distance between each atom in the interface  
        energy=0
        
        m1=self.multimer.get_multimer_uxyz()[0]
        m2=self.multimer.get_multimer_uxyz()[1]
        
        # list of the residue name, atomname, chain and residue number
        data_list = np.array(self.data.structure.mapping(self.data.structure.data))[:,1:5]
              
        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
        dist=np.array(d)

        
        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #the indexes of the clashing atoms are stored here
        
        
        # extract all the atoms of the amino acids involved in the interface:
        interface1 = []
        interface2 = []
        
        # due to the fact that sometimes pdb residues have different nomenclatures than those of CHARMM27 rtf
        # we need to specify a dictionary enabling the conversion
        resDict = {"HIS":"HSE"
                   }
        
        # we want to take all the atoms contained the residue where at least one atom is presen the interface (couples) index
        for i in xrange(0,len(couples[0]),1):
            
            if couples[0][i] not in interface1:
                # extract the interface atom RESIDUE number
                residueNumber1=self.data.data_list[couples[0][i]][-1] 
                # add all the ATOMS contained in the RESIDUE containing at least one atom involved in the first interface                            
                interface1 = interface1 + list(np.where(self.data.data_list[:,-1] == residueNumber1)[0])
                                
            if couples[1][i] not in interface2:
                residueNumber2=self.data.data_list[couples[1][i]][-1]   
                interface2 = interface2 + list(np.where(self.data.data_list[:,-1] == residueNumber2)[0])
   
        for atom1 in interface1:
            for atom2 in interface2:       
                # storing the information about each atom, namely the residue they are contained in and the atom assigned to them
                atom1Resname = self.data.data_list[atom1][1]
                atom2Resname = self.data.data_list[atom2][1]
                atom1PDBname = self.data.data_list[atom1][0]
                atom2PDBname = self.data.data_list[atom2][0]
                
                
                # calculating VDW iteractions
                rij = dist[atom1,atom2]
                
                Ei = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][0]][0]
                Ej = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][0]][0]
                Eij = np.sqrt( Ei * Ej )
                
                ri = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][0]][1]
                rj = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][0]][1]
                rmin = ri+rj
                
                VW = Eij*( (rmin/rij)**12 - (2*(rmin/rij)**6) )
                
                # calculating electrostatic interactions
                K = 332.0637
                qi = self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][1]
                qj = self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][1]
                
                Vele = K * ( (qi * qj) / rij )
                
                # adding all to calculate the non bonded energy
                energy += Vele + VW
                    
        return energy
    
    def AllAtomNonVW(self):
        cutoff = 3.0 # the distance between each atom in the interface  
        energy=0
        
        m1=self.multimer.get_multimer_uxyz()[0]
        m2=self.multimer.get_multimer_uxyz()[1]
                      
        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
        dist=np.array(d)

        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #the indexes of the clashing atoms are stored here
        
        
        # extract all the atoms of the amino acids involved in the interface:
#         interface1 = []
#         interface2 = []
        
        # we want to take all the atoms contained the residue where at least one atom is presen the interface (couples) index
#         for i in xrange(0,len(couples[0]),1):
#             
#             if couples[0][i] not in interface1:
#                 # extract the interface atom RESIDUE number
#                 residueNumber1=self.data.data_list[couples[0][i]][-1] 
#                 # add all the ATOMS contained in the RESIDUE containing at least one atom involved in the first interface                            
#                 interface1 = interface1 + list(np.where(self.data.data_list[:,-1] == residueNumber1)[0])
#                                 
#             if couples[1][i] not in interface2:
#                 residueNumber2=self.data.data_list[couples[1][i]][-1]   
#                 interface2 = interface2 + list(np.where(self.data.data_list[:,-1] == residueNumber2)[0])
   
        for i in xrange(0,len(couples[0]),1):
            atom1 = couples[0][i]
            atom2 = couples[1][i]
            
            # storing the information about each atom, namely the residue they are contained in and the atom assigned to them
            atom1Resname = self.data.data_list[atom1][1]
            atom2Resname = self.data.data_list[atom2][1]
            atom1PDBname = self.data.data_list[atom1][0]
            atom2PDBname = self.data.data_list[atom2][0]
            
            # calculating VDW iteractions
            rij = dist[atom1,atom2]
            
            Ei = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][0]][0]
            Ej = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][0]][0]
            Eij = np.sqrt( Ei * Ej )
            
            ri = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][0]][1]
            rj = self.data.AtomVDWprm[ self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][0]][1]
            rmin = ri+rj
            
            VW = Eij*( (rmin/rij)**9 - (2*(rmin/rij)**6) )
            
            # calculating electrostatic interactions
#             K = 332.0637
#             qi = self.data.ResidueAtomCharge[atom1Resname][atom1PDBname][1]
#             qj = self.data.ResidueAtomCharge[atom2Resname][atom2PDBname][1]
#             
#             Vele = K * ( (qi * qj) / rij )
            
            # adding all to calculate the non bonded energy
#             energy += Vele + VW
            energy += VW
        
        return energy
    
    def CGNonBonded(self,pos):
# first create the CG model of the multimer:
#         pdb_file = self.multimer.internal_PDB()
#         
#         cgData = martini.runFromPOW(['-f', pdb_file, '-x', 'null_cg.pdb', '-p', 'backbone', '-ff', 'martini22'])
#         self.cg_multimer = M.Multimer(self.data.cg_structure)
#         self.cg_multimer.create_multimer(self.degree,pos[3],pos[0:3])
         
        m1=self.multimer.get_multimer_cg_uxyz()[0]
        m2=self.multimer.get_multimer_cg_uxyz()[1]
        
#         m1 = np.array(cgData[0])
#         m2 = np.array(cgData[1])
 
        epsilon=1.0
        sigma=4.7
        cutoff=12.0

        energy=0

        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

        dist=np.array(d)

        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues

        for i in xrange(0,len(couples[0]),1):
            d=dist[couples[0,i],couples[1,i]]
            energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6) ## was 9
        
        return energy
        
        

class Postprocess(PP):

    def __init__(self,data,params):
        self.data=data
        self.params=params

        self.rec_dim=0
        if params.receptor!="NA":
            self.rec_dim=4

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
        
        import ClusterAndDraw as CnD
        import wx
        
        
        multimer = M.Multimer(self.data.structure)
        multimer.create_multimer(2,10,np.array([0.0,0.0,0.0]))
        [m,self.index]=multimer.atomselect(1,"*","*","CA",True) # -> extracting indexes of CA

        #load the monomeric structure positions (needed for resetting atom position after displacement)
        s = Protein()
        s.import_pdb(self.params.pdb_file_name)
        coords=s.get_xyz()
        self.coords = deepcopy(coords)
        
        Matrix_creator = CnD.Matrix_creator(self.params, self.data, self)
        Matrix_creator.computeMatrix()
        
        comm.Barrier()
        self.OUTPUT_DIRECTORY=self.params.output_folder
        if rank == 0:
            
            if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
                os.mkdir(self.OUTPUT_DIRECTORY)
            
            #self.OUTPUT_DIRECTORY=self.params.output_folder
            app = wx.App(False)
            frame = CnD.MainFrame(None, "Clustering interface",self.OUTPUT_DIRECTORY ,self.params, self.data, self)
            frame.distancePanel.computeCluster()
            
            if self.params.cluster_threshold == "NA":
                frame.Show()
                app.MainLoop()
            else:
                frame.distancePanel.convertCoordsAndExport(self.params.cluster_threshold)

    #this function superseed the standard euclidean distance implemented in Default.py
    def computeDistance (self, data1, data2):            
        
        # create the first multimer
        if self.params.style=="flexible":
            #pick monomeric structure from database
            deform_coeffs=data1[(4+self.rec_dim):-1]
            #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
            #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
            code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
            target_frame=min_dist.argmin()
            coords=self.data.traj[:,target_frame]
            coords_reshaped=coords.reshape(len(coords)/3,3)
            self.data.structure.set_xyz(coords_reshaped)                         
        else:
            self.data.structure.set_xyz(self.coords)

        #pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
        multimer1 = M.Multimer(self.data.structure)
        multimer1.create_multimer(self.params.degree,data1[3],np.array([data1[0],data1[1],data1[2]]))
        m1=multimer1.get_multimer_uxyz()[0][self.index] # getting all the coordinates of m1

        #create the second multimer
        if self.params.style=="flexible":
            #pick monomeric structure from database
            deform_coeffs=data2[(4+self.rec_dim):-1]
            #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
            #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
            code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
            target_frame=min_dist.argmin()
            coords=self.data.traj[:,target_frame]
            coords_reshaped=coords.reshape(len(coords)/3,3)
            self.data.structure.set_xyz(coords_reshaped)                         
        else:
            self.data.structure.set_xyz(self.coords)
            
        multimer2 = M.Multimer(self.data.structure)
        multimer2.create_multimer(self.params.degree,data2[3],np.array([data2[0],data2[1],data2[2]]))
        m2=multimer2.get_multimer_uxyz()[0][self.index]

        # calculate distance between the 2
        rmsd=self.align(m1,m2) # --> comes from Default.Postprocess.align()
        
        return rmsd
                
                
    def make_output(self, coordinateArray, average_RMSD_ARRAY):
        
        self.coordinateArray = coordinateArray
        print self.coordinateArray 
                    
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

        s = Protein()
        s.import_pdb(self.params.pdb_file_name)
        coords=s.get_xyz()
        
        for n in xrange(0,coordinateArray.shape[0],1):
            
            print "creating PDB for structure index: "+str(n)
            # create the first multimer
            if self.params.style=="flexible":
                #pick monomeric structure from database
                deform_coeffs=coordinateArray[n][(4+self.rec_dim):-1]
                #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
                #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
                code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
                target_frame=min_dist.argmin()
                coords=self.data.traj[:,target_frame]
                coords_reshaped=coords.reshape(len(coords)/3,3)
                self.data.structure.set_xyz(coords_reshaped)                         
            else:
                self.data.structure.set_xyz(coords)
    
            #pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            multimer1 = M.Multimer(self.data.structure)
            multimer1.create_multimer(self.params.degree, coordinateArray[n][3], np.array([coordinateArray[n][0],coordinateArray[n][1],coordinateArray[n][2]]))

            # print the pdb file
            multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,iterant+1))

            # create the constraint:
            measure = constraint.constraint_check(multimer1)

            # WRITING SOLUTION.DAT
            l = []
            f = []

            # insert coordinates in the solution.dat file
            f.append("assembly "+str(iterant+1)+" |")
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
            # check here if the structure belongs to a cluster
            l.append(average_RMSD_ARRAY[iterant])

            formatting=''.join(f)

            clusters_file.write(formatting%tuple(l))

            # --------------------------- WRITING TCL FILE
            if iterant == 0:
                tcl_file.write("mol new assembly"+str(iterant+1)+".pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all \n")

            else:
                tcl_file.write("mol addfile assembly"+str(iterant+1)+".pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all \n")

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
        
        
    def run2(self) :
        
        exec 'import %s as constraint'%(self.constraint)

        #create output directory for generated PDB
        self.OUTPUT_DIRECTORY=self.params.output_folder
        if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
            os.mkdir(self.OUTPUT_DIRECTORY)

        clusters_file=open("%s/solutions.dat"%self.params.output_folder,"w")

        #use superclass method to filter acceptable solutions
        self.log=self.select_solutions(self.params)
        print ">> %s solutions filtered"%len(self.log[:,1])
        if len(self.log[:,1])==0:
            return
    
        #generate a dummy multimer and extract the indexes of C alpha
        multimer = M.Multimer(self.data.structure)
        multimer.create_multimer(2,10,np.array([0.0,0.0,0.0]))
        [m,index]=multimer.atomselect(1,"*","*","CA",True)

        #if needed, extract receptor indexes too
        if self.data.receptor_structure!=[]:
            [r1,index_receptor]=self.data.receptor_structure.atomselect("*","*","CA",True)

        #load the monomeric structure positions (needed for resetting atom position after displacement)
        s = Protein()
        s.import_pdb(self.params.pdb_file_name)   
        coords=s.get_xyz()
        self.coords = deepcopy(coords)

        print ">> clustering best solutions..."
        P=self.log[:,0:len(self.log[0,:])] #points list
        V=self.log[:,-1] #values of best hits
        C=[] #centroids array
        P_C=np.zeros(len(P)) #points-to-cluster mapping
        C_V=[] #centroids values
        cnt=0 #centroids counter

        #cluster accepted solutions
        while(True) :
            #check if new clustering loop is needed
            k=np.nonzero(P_C==0)[0]
            if len(k)!=0 :
                cnt=cnt+1
                P_C[k[0]]=cnt
                a=P[k[0]]
                C.append(a)
            else :
                break

            if self.params.style=="flexible":
                #build reference structure
                deform_coeffs=np.array(C[cnt-1])[(4+self.rec_dim):len(C[cnt-1])-1]

                pos=self.data.proj[:,self.data.centroid]+deform_coeffs
                code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
                target_frame1=min_dist.argmin()
                coords=self.data.traj[:,target_frame1]
                coords_reformat=coords.reshape(len(coords)/3,3)
                self.data.structure.set_xyz(coords_reformat)
            else:
                self.data.structure.set_xyz(coords)  

            #create multimer
            pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            print pos
            multimer1 = M.Multimer(self.data.structure)
            multimer1.create_multimer(self.params.degree,pos[3],pos[0:3])
            m1=multimer1.get_multimer_uxyz()[0][index]

            #if self.data.receptor_structure!=[]:
                ##shift multimer on z axis
                #multimer1.z_to_origin()
                #multimer1.z_shift(pos[4])

                ##rotate receptor (backing up original position, and putting it back when measure is done),
                ##also write the rotated receptor structure
                #coords_tmp=self.data.receptor_structure.get_xyz()
                #self.data.receptor_structure.rotation(pos[5],pos[6],pos[7])
                #r1=self.data.receptor_structure.get_xyz()[index_receptor]
                #self.data.receptor_structure.write_pdb("%s/receptor%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))
                #self.data.receptor_structure.set_xyz(coords_tmp)
                #m1=np.concatenate((m1,r1))

            #write multimer
            multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))

            #clustering loop
            cnt2=1
            for i in xrange(0,len(k),1) :

                if self.params.style=="flexible":
                    deform_coeffs=np.array(P[k[i]])[(4+self.rec_dim):len(P[k[i]])-1]
                    pos=self.data.proj[:,self.data.centroid]+deform_coeffs
                    code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
                    target_frame=min_dist.argmin()
                    coords=self.data.traj[:,target_frame]
                    coords_reformat=coords.reshape(len(coords)/3,3)
                    self.data.structure.set_xyz(coords_reformat)
                else:
                    self.data.structure.set_xyz(coords)  

                multimer2 = M.Multimer(self.data.structure)
                multimer2.create_multimer(self.params.degree,P[k[i]][3],np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2]]))
                m2=multimer2.get_multimer_uxyz()[0][index]

                if self.data.receptor_structure!=[]:
                    #shift multimer on z axis
                    multimer2.z_to_origin()
                    multimer2.z_shift(P[k[i]][4])

                    #rotate receptor (backing up original position, and putting it back when measure is done)
                    coords_tmp=self.data.receptor_structure.get_xyz()
                    self.data.receptor_structure.rotation(P[k[i]][5],P[k[i]][6],P[k[i]][7])
                    r2=self.data.receptor_structure.get_xyz()[index_receptor]
                    self.data.receptor_structure.set_xyz(coords_tmp)
                    m2=np.concatenate((m2,r2))

                #compute RMSD
                rmsd=self.align(m1,m2)
                            
                if rmsd<self.params.cluster_threshold :
                    cnt2+=1
                    P_C[k[i]]=cnt
            
            if self.params.style=="rigid":
                print ">>> clustered %s solutions on multimer %s"%(cnt2-1,cnt)
            if self.params.style=="flexible":
                print ">>> clustered %s solutions on multimer %s (from frame %s)"%(cnt2-1,cnt,target_frame1)
      
            #set centroid score with score of closes neighbor in set
            q=np.nonzero(P_C==cnt)[0]
            distance=10000
            targ=0
            for i in xrange(0,len(q),1) :
                d=np.sqrt(np.dot(C[cnt-1]-P[q[i]],C[cnt-1]-P[q[i]]))         
                if d<distance :
                    distance=d
                    targ=q[i]
            C_V.append(V[targ])
            
            #extract constraint values calculated for selected centroid
            measure = constraint.constraint_check(multimer1)

            ###generate output log (prepare data and formatting line, then dump in output file)###
            l=[]
            f=[]
            for item in C[cnt-1][0:len(C[cnt-1])-1]:
                l.append(item)
                f.append("%8.3f ")
            #write constraint values
            f.append("| ")
            for item in measure:
                l.append(item)
                f.append("%8.3f ")
            #write fitness
            f.append("| %8.3f\n")
            l.append(C_V[cnt-1])

            formatting=''.join(f)

            clusters_file.write(formatting%tuple(l))


        clusters_file.close()

        
    
    def run3(self):
                
        if rank == 0:
        
            #create output directory for generated PDB
            self.OUTPUT_DIRECTORY=self.params.output_folder
            if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
                os.mkdir(self.OUTPUT_DIRECTORY)
              
            iterant = 0 # this iterator is used count the number of assemblies created
            
            # open the whole log file to parse through
            logFile = open(self.params.output_file)
            temp = []
            
            #read logfile into memory and extract coordinates only:
            while True:
                line = logFile.readline()
                if line:
                    temp.append(line.split(" "))
                else:
                    break
                
           
                
            coordinateArray = np.array(temp)[:,2:6].astype("float")
            
    
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
    
            s = Protein()
            s.import_pdb(self.params.pdb_file_name)
            coords=s.get_xyz()
            
            
        else:
            s = None
            coords=None
            coordinateArray=None
            self.OUTPUT_DIRECTORY=None
            
        comm.Barrier()
        s=comm.bcast( s,root=0)
        coords=comm.bcast( coords,root=0)
        coordinateArray=comm.bcast(coordinateArray, root=0)
        self.OUTPUT_DIRECTORY=comm.bcast(self.OUTPUT_DIRECTORY, root=0)
        comm.Barrier()
            
        n=rank
        
        
#        for n in xrange(0,coordinateArray.shape[0],1):
        while n < coordinateArray.shape[0]:    
            # create the first multimer
            if self.params.style=="flexible":
                #pick monomeric structure from database
                deform_coeffs=coordinateArray[n][(4+self.rec_dim):-1]
                #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
                #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
                code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
                target_frame=min_dist.argmin()
                coords=self.data.traj[:,target_frame]
                coords_reshaped=coords.reshape(len(coords)/3,3)
                self.data.structure.set_xyz(coords_reshaped)                         
            else:
                self.data.structure.set_xyz(coords)
    
            #pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            multimer1 = M.Multimer(self.data.structure)
            multimer1.create_multimer(self.params.degree, coordinateArray[n][3], np.array([coordinateArray[n][0],coordinateArray[n][1],coordinateArray[n][2]]))

            # print the pdb file
            if 11 < 3:
                multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,iterant+1))
#                iterant += 1
            else:
                multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,n+1))
                n += size
                
            if rank == 0:
                iterant += 1
                
                if n % 100 == 0:
                    print "creating PDB for structure index: "+str(n)
            
        print "process "+str(rank)+" finished normally."
        comm.Barrier()
               
        if rank == 0:   
            logFile.close()  
            print ">> Done make .pdbs"            
            
            
    def run_cg(self) :
        
        exec 'import %s as constraint'%(self.constraint)

        #create output directory for generated PDB
        self.OUTPUT_DIRECTORY=self.params.output_folder
        if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
            os.mkdir(self.OUTPUT_DIRECTORY)

        clusters_file=open("%s/solutions.dat"%self.params.output_folder,"w")

        #use superclass method to filter acceptable solutions
        self.log=self.select_solutions(self.params)
        print ">> %s solutions filtered"%len(self.log[:,1])
        if len(self.log[:,1])==0:
            return
    
        #generate a dummy multimer and extract the indexes of C alpha
        multimer = MCG.Multimer(self.data.structure, self.data.cg_structure)
        multimer.create_multimer(2,10,np.array([0.0,0.0,0.0]))
        [m,index]=multimer.atomselect(1,"*","*","CA",True)
        
        

        #if needed, extract receptor indexes too
        if self.data.receptor_structure!=[]:
            [r1,index_receptor]=self.data.receptor_structure.atomselect("*","*","CA",True)

        #load the monomeric structure positions (needed for resetting atom position after displacement)
        s = Protein()
        s.import_pdb(self.params.pdb_file_name)   
        coords=s.get_xyz()
        
        cg = Protein()
        cg.import_pdb('CG_protein_0.pdb')
        cg_coords = cg.get_xyz()
        

        print ">> clustering best solutions..."
        P=self.log[:,0:len(self.log[0,:])] #points list
        V=self.log[:,-1] #values of best hits
        C=[] #centroids array
        P_C=np.zeros(len(P)) #points-to-cluster mapping
        C_V=[] #centroids values
        cnt=0 #centroids counter

        #cluster accepted solutions
        while(True) :
            #check if new clustering loop is needed
            k=np.nonzero(P_C==0)[0]
            if len(k)!=0 :
                cnt=cnt+1
                P_C[k[0]]=cnt
                a=P[k[0]]
                C.append(a)
            else :
                break

#             if self.params.style=="flexible":
#                 #build reference structure
#                 deform_coeffs=np.array(C[cnt-1])[(4+self.rec_dim):len(C[cnt-1])-1]
# 
#                 pos=self.data.proj[:,self.data.centroid]+deform_coeffs
#                 code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
#                 target_frame1=min_dist.argmin()
#                 coords=self.data.traj[:,target_frame1]
#                 coords_reformat=coords.reshape(len(coords)/3,3)
#                 self.data.structure.set_xyz(coords_reformat)
#             else:
            self.data.structure.set_xyz(coords)
            self.data.cg_structure.set_xyz(cg_coords)

            #create multimer
            pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            multimer1 = MCG.Multimer(self.data.structure, self.data.cg_structure)
            multimer1.create_multimer(self.params.degree,pos[3],pos[0:3])
            m1=multimer1.get_multimer_uxyz()[0][index]
            
#             multimer_cg = M.Multimer(self.data.cg_structure)
#             multimer_cg.create_multimer(self.params.degree,pos[3],pos[0:3])

            #if self.data.receptor_structure!=[]:
                ##shift multimer on z axis
                #multimer1.z_to_origin()
                #multimer1.z_shift(pos[4])

                ##rotate receptor (backing up original position, and putting it back when measure is done),
                ##also write the rotated receptor structure
                #coords_tmp=self.data.receptor_structure.get_xyz()
                #self.data.receptor_structure.rotation(pos[5],pos[6],pos[7])
                #r1=self.data.receptor_structure.get_xyz()[index_receptor]
                #self.data.receptor_structure.write_pdb("%s/receptor%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))
                #self.data.receptor_structure.set_xyz(coords_tmp)
                #m1=np.concatenate((m1,r1))

            #write multimer
            multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))
            multimer1.write_CG("%s/cg_assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))

            #clustering loop
            cnt2=1
            for i in xrange(0,len(k),1) :

                if self.params.style=="flexible":
                    deform_coeffs=np.array(P[k[i]])[(4+self.rec_dim):len(P[k[i]])-1]
                    pos=self.data.proj[:,self.data.centroid]+deform_coeffs
                    code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
                    target_frame=min_dist.argmin()
                    coords=self.data.traj[:,target_frame]
                    coords_reformat=coords.reshape(len(coords)/3,3)
                    self.data.structure.set_xyz(coords_reformat)
                else:
                    self.data.structure.set_xyz(coords)  

                multimer2 = M.Multimer(self.data.structure)
                multimer2.create_multimer(self.params.degree,P[k[i]][3],np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2]]))
                m2=multimer2.get_multimer_uxyz()[0][index]

                if self.data.receptor_structure!=[]:
                    #shift multimer on z axis
                    multimer2.z_to_origin()
                    multimer2.z_shift(P[k[i]][4])

                    #rotate receptor (backing up original position, and putting it back when measure is done)
                    coords_tmp=self.data.receptor_structure.get_xyz()
                    self.data.receptor_structure.rotation(P[k[i]][5],P[k[i]][6],P[k[i]][7])
                    r2=self.data.receptor_structure.get_xyz()[index_receptor]
                    self.data.receptor_structure.set_xyz(coords_tmp)
                    m2=np.concatenate((m2,r2))

                #compute RMSD
                rmsd=self.align(m1,m2)
                            
                if rmsd<self.params.cluster_threshold :
                    cnt2+=1
                    P_C[k[i]]=cnt
            
            if self.params.style=="rigid":
                print ">>> clustered %s solutions on multimer %s"%(cnt2-1,cnt)
            if self.params.style=="flexible":
                print ">>> clustered %s solutions on multimer %s (from frame %s)"%(cnt2-1,cnt,target_frame1)
      
            #set centroid score with score of closes neighbor in set
            q=np.nonzero(P_C==cnt)[0]
            distance=10000
            targ=0
            for i in xrange(0,len(q),1) :
                d=np.sqrt(np.dot(C[cnt-1]-P[q[i]],C[cnt-1]-P[q[i]]))         
                if d<distance :
                    distance=d
                    targ=q[i]
            C_V.append(V[targ])
            
            #extract constraint values calculated for selected centroid
            measure = constraint.constraint_check(multimer1)

            ###generate output log (prepare data and formatting line, then dump in output file)###
            l=[]
            f=[]
            for item in C[cnt-1][0:len(C[cnt-1])-1]:
                l.append(item)
                f.append("%8.3f ")
            #write constraint values
            f.append("| ")
            for item in measure:
                l.append(item)
                f.append("%8.3f ")
            #write fitness
            f.append("| %8.3f\n")
            l.append(C_V[cnt-1])

            formatting=''.join(f)

            clusters_file.write(formatting%tuple(l))


        clusters_file.close()
        
        
    def run_clusco(self):
                
        if rank == 0:
        
            #create output directory for generated PDB
            self.OUTPUT_DIRECTORY=self.params.output_folder
            if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
                os.mkdir(self.OUTPUT_DIRECTORY)
              
            iterant = 0 # this iterator is used count the number of assemblies created
            
            # open the whole log file to parse through
            logFile = open(self.params.output_file)
            temp = []
            
            #read logfile into memory and extract coordinates only:
            while True:
                line = logFile.readline()
                if line:
                    temp.append(line.split(" "))
                else:
                    break
                
           
                
            coordinateArray = np.array(temp)[:,2:6].astype("float")
            
    
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
    
            s = Protein()
            s.import_pdb(self.params.pdb_file_name)
            coords=s.get_xyz()
            
            
        else:
            s = None
            coords=None
            coordinateArray=None
            self.OUTPUT_DIRECTORY=None
            
        comm.Barrier()
        s=comm.bcast( s,root=0)
        coords=comm.bcast( coords,root=0)
        coordinateArray=comm.bcast(coordinateArray, root=0)
        self.OUTPUT_DIRECTORY=comm.bcast(self.OUTPUT_DIRECTORY, root=0)
        comm.Barrier()
            
        n=rank
        
        
#        for n in xrange(0,coordinateArray.shape[0],1):
        while n < coordinateArray.shape[0]:    
            # create the first multimer
            if self.params.style=="flexible":
                #pick monomeric structure from database
                deform_coeffs=coordinateArray[n][(4+self.rec_dim):-1]
                #pos_eig=self.data.proj[:,self.data.centroid]+deform_coeffs
                #code,min_dist=vq(self.data.proj.transpose(),np.array([pos_eig]))
                code,min_dist=vq(self.data.proj.transpose(),np.array([deform_coeffs]))
                target_frame=min_dist.argmin()
                coords=self.data.traj[:,target_frame]
                coords_reshaped=coords.reshape(len(coords)/3,3)
                self.data.structure.set_xyz(coords_reshaped)                         
            else:
                self.data.structure.set_xyz(coords)
    
            #pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            multimer1 = M.Multimer(self.data.structure)
            multimer1.create_multimer(self.params.degree, coordinateArray[n][3], np.array([coordinateArray[n][0],coordinateArray[n][1],coordinateArray[n][2]]))

            # print the pdb file
            multimer1.write_PDB_clusco("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,n+1))
            n += size
                
            if rank == 0:
                iterant += 1
                
                if n % 100 == 0:
                    print "creating PDB for structure index: "+str(n)
            
        print "process "+str(rank)+" finished normally."
        comm.Barrier()
               
        if rank == 0:   
            logFile.close()  
            print ">> Done make .pdbs" 
            
            
        ###########################################################################
    
    
    
    def select_solutions(self, logArray, logFileType):
            '''- returns indeces of converged solution in the case of mVie
               - returns paramters of top n % assemblies in case of pow/pso
            '''
            
            selected = []
            goodIndeces = []
            flag = 1
            
            # case of mvie
            if logFileType == 0:
                # keep indeces where boundaries are converged
                for i in xrange(0,len(logArray,),1):
                    flag = 1
                    for j in xrange(0,len(logArray[i][5:]),1):
                        if logArray[i][5:][j] != -2:
                            flag = 0
                            break
                    if flag == 1:
                        goodIndeces.append(i)
                return goodIndeces
            # case of pso/pow
            else:
                # sort logArray according to fitness and get indeces of top n%
                a = sorted(logArray, key=lambda a_entry: a_entry[-1]) 
                b = np.array(a)
                
                topN = int((25.0 * float(len(b)))/100.0)
                top  = b[0:topN]
                
                return top
                
             
            
       
    def get_centroids(self, logFileName, logFileType, destFolder, cccFlag) :
        '''Unlike run_2(), takes into account the rmsd associated with 
           assemblies. Also differentiate between mVie and pso file types.
        '''
        
        #create output directory for generated PDB
        self.OUTPUT_DIRECTORY = destFolder
        if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
            os.mkdir(self.OUTPUT_DIRECTORY)
          
        iterant = 0 # this iterator is used count the number of assemblies created
        
        # open the whole log file to parse through
        logFile = open(logFileName)
        temp = []
        
        if logFileType == 0 :
           separator="\t"
        else :
           separator=" "

        # --- read logfile into memory and extract coordinates only:
        while True:
            line = logFile.readline()
            if line:
                temp.append(line.split(separator))
            else:
                break
        
        # --- extract the select assemblies depending on type of algo used
        if logFileType == 0 : # -> mVie
            # here we use the function above to extract the indeces of 
            # solutions satisfying spatial restraints
            selected = self.select_solutions(np.array(temp)[:,0:5+(len(self.params.target)*2)].astype("float"),logFileType)
            self.log=np.array(temp)[:,0:4].astype("float")[selected]
        
        else : # -> PSO
            # select n% best assemblies based on fitness function
            selected = self.select_solutions(np.array(temp)[:,2:].astype("float"),logFileType)
            self.log = selected[:,0:-1]
            
            # extract rmsds from rmsd files and add values to temp
            rmsdFile = open(logFileName[:-13]+"_rmsd.txt")
            
            # add the rmsds to the coordinates of the assemblies
            # this step was not necessary above as the logfile of mVie already
            # contains rmsds associated with parameters
            npTemp = np.array(temp)
            rmsds = np.array([[line[:-1] for line in rmsdFile]])
            temp = np.concatenate((npTemp, rmsds.T), axis=1)[:,2:]
            
            
            
                    
        print ">> %s solutions filtered"%len(self.log[:,1])
        if len(self.log[:,1])==0:
            return
        
        # create dictionary having for key the parameters and value rmsds
        # -> necessary to assign rmsds to centroids
        paramsToRmsd = {}
        if logFileType == 0 :
            for item in np.array(temp).astype("float")[selected]:
                paramsToRmsd[" ".join([str(x) for x in item[0:4]])] = str(item[4])+" "+str(item[-1])
#                 paramsToRmsd[" ".join([str(x) for x in item[0:4]])] = str(item[-1])
        else:
            for item in np.array(temp).astype("float"):
                paramsToRmsd[" ".join([str(x) for x in item[0:4]])] = " ".join( [str(x) for x in item[4:]] )
#                 paramsToRmsd[" ".join([str(x) for x in item[0:4]])] = str(item[-1])
            
        # open file for writing the centroids parameters and associated rmsd
        clusterFile = open(self.OUTPUT_DIRECTORY+"/"+logFileName.split("/")[-1],"w")
        
        
        #generate a dummy multimer and extract the indexes of C alpha
        multimer = M.Multimer(self.data.structure)
        multimer.create_multimer(2,10,np.array([0.0,0.0,0.0]))
        [m,index]=multimer.atomselect(1,"*","*","CA",True)

        #if needed, extract receptor indexes too
        if self.data.receptor_structure!=[]:
            [r1,index_receptor]=self.data.receptor_structure.atomselect("*","*","CA",True)

        #load the monomeric structure positions (needed for resetting atom position after displacement)
        s = Protein()
        s.import_pdb(self.params.pdb_file_name)   
        coords=s.get_xyz()
        self.coords = deepcopy(coords)

        print ">> clustering best solutions..."
        P=self.log[:,0:len(self.log[0,:])] #points list
        V=self.log[:,-1] #values of best hits
        C=[] #centroids array
        P_C=np.zeros(len(P)) #points-to-cluster mapping
        C_V=[] #centroids values
        cnt=0 #centroids counter

        #cluster accepted solutions
        while(True) :
            #check if new clustering loop is needed
            k=np.nonzero(P_C==0)[0]
            if len(k)!=0 :
                cnt=cnt+1
                P_C[k[0]]=cnt
                a=P[k[0]]
                C.append(a)
            else :
                break

            if self.params.style=="flexible":
                #build reference structure
                deform_coeffs=np.array(C[cnt-1])[(4+self.rec_dim):len(C[cnt-1])-1]

                pos=self.data.proj[:,self.data.centroid]+deform_coeffs
                code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
                target_frame1=min_dist.argmin()
                coords=self.data.traj[:,target_frame1]
                coords_reformat=coords.reshape(len(coords)/3,3)
                self.data.structure.set_xyz(coords_reformat)
            else:
                self.data.structure.set_xyz(coords)  

            #create multimer
            pos = np.array(C[cnt-1])[0:(4+self.rec_dim)].astype(float)
            
            # from position of cluster, extract rmsd:
            key = " ".join([str(x) for x in pos])
            print str(key)+" >>> "+str(paramsToRmsd[key])
            clusterFile.write(str(key)+" "+str(paramsToRmsd[key]) )
            
            
            multimer1 = M.Multimer(self.data.structure)
            multimer1.create_multimer(self.params.degree,pos[3],pos[0:3])
            m1=multimer1.get_multimer_uxyz()[0][index]

            # ccc extraction if flag is on:
            if cccFlag == 1:
                multimer1.create_PDB_for_density_map(rank, self.params.pdb_file_name)
                ccc = 1 - CCC.get_DMAP_fit( str(rank), self.params.map_resolution ) 
                print "> 1 - ccc: "+str(ccc)
                clusterFile.write(" "+str(ccc))
                
            clusterFile.write("\n")
            
            #write multimer
#             multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))
            
            #clustering loop
            cnt2=1
            for i in xrange(0,len(k),1) :

                if self.params.style=="flexible":
                    deform_coeffs=np.array(P[k[i]])[(4+self.rec_dim):len(P[k[i]])-1]
                    pos=self.data.proj[:,self.data.centroid]+deform_coeffs
                    code,min_dist=vq(self.data.proj.transpose(),np.array([pos]))
                    target_frame=min_dist.argmin()
                    coords=self.data.traj[:,target_frame]
                    coords_reformat=coords.reshape(len(coords)/3,3)
                    self.data.structure.set_xyz(coords_reformat)
                else:
                    self.data.structure.set_xyz(coords)  

                multimer2 = M.Multimer(self.data.structure)
                multimer2.create_multimer(self.params.degree,P[k[i]][3],np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2]]))
                m2=multimer2.get_multimer_uxyz()[0][index]

                #compute RMSD
                rmsd=self.align(m1,m2)
                            
                if rmsd<self.params.cluster_threshold :
                    cnt2+=1
                    P_C[k[i]]=cnt
            
#             if self.params.style=="rigid":
#                 print ">>> clustered %s solutions on multimer %s"%(cnt2-1,cnt)
            if self.params.style=="flexible":
                print ">>> clustered %s solutions on multimer %s (from frame %s)"%(cnt2-1,cnt,target_frame1)
      
            #set centroid score with score of closes neighbor in set
            q=np.nonzero(P_C==cnt)[0]
            distance=10000
            targ=0
            for i in xrange(0,len(q),1) :
                d=np.sqrt(np.dot(C[cnt-1]-P[q[i]],C[cnt-1]-P[q[i]]))         
                if d<distance :
                    distance=d
                    targ=q[i]
            C_V.append(V[targ])
            
        
        clusterFile.close()

            
        
        
