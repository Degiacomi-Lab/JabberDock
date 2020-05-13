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
# Author : Matteo Degiacomi, matteothomas.degiacomi@epfl.ch# Copyright (c) 2012 EPFL (Ecole Polytechnique federale de Lausanne)
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
import Assembly as A
import flexibility_new as F


class Parser(R):

    def __init__(self):

        self.add('rmsdSelect','rmsd_select','array str',"NA")
        self.add('constraint','constraint','str',"NA")
        self.add('energy','energy_type','str',"vdw")
        self.add('detectClash','detect_clash','str',"on")
        self.add('target','target','array float',np.array([]))
        self.add('mode','mode','str',"seed")

        #ligand flags
        self.add('ligand_style','ligand_style','str',"rigid")
        self.add('ligand_projection','ligand_proj_file','str',"NA")
        self.add('ligand_align','ligand_align','str',"no")
        self.add('ligand_ratio','ligand_ratio','float',0.9)
        self.add('ligand_trajectory','ligand_trajectory','str',"NA")
        self.add('ligand_trajSelection','ligand_trajselection','str',"NA")
        self.add('ligand_topology','ligand_topology','str',"NA")
        self.add('ligand','ligand_file_name','str',"NA")

        #receptor flags
        self.add('receptor_style','receptor_style','str',"rigid")
        self.add('receptor_projection','receptor_proj_file','str',"NA")
        self.add('receptor_align','receptor_align','str',"no")
        self.add('receptor_ratio','receptor_ratio','float',0.9)
        self.add('receptor_trajectory','receptor_trajectory','str',"NA")
        self.add('receptor_trajSelection','receptor_trajselection','str',"NA")
        self.add('receptor_topology','receptor_topology','str',"NA")
        self.add('receptor','receptor_file_name','str',"NA")

        self.add('cluster_threshold','cluster_threshold','float',"2")
        self.add('output_folder','output_folder','str',"result")


    def check_variables(self):

        if self.cluster_threshold<0:
            print "ERROR: clustering threshlod should be greater than 0!"
            sys.exit(1)

        #check ligand files existence
        if self.ligand_file_name!="NA" and self.ligand_style=="rigid":
            tmp=os.path.abspath(self.ligand_file_name)
            if os.path.isfile(self.ligand_file_name)!=1 :
                print "ERROR: ligand pdb file %s not found"%self.ligand_file_name
                sys.exit(1)

        #check receptor existence
        if self.receptor_file_name!="NA" and self.receptor_style=="rigid":
            tmp=os.path.abspath(self.receptor_file_name)
            if os.path.isfile(self.receptor_file_name)!=1 :
                print "ERROR: receptor pdb file %s not found"%self.receptor_file_name
                sys.exit(1)

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



class Data:

    index_ligand=[]
    index_receptor=[]
    cg_atoms=[]

    def __init__(self,params):

        #LIGAND STRUCTURE
        self.ligand = Protein()

        if params.ligand_style=="flexible":
            print ">> flexible docking requested for ligand, launching PCA..."
            try:
                self.flex_ligand=F.Flexibility_PCA()
                self.flex_ligand.compute_eigenvectors(params.ligand_topology,params.ligand_trajectory,params.ligand_align,params.ligand_ratio,params.mode,params.ligand_proj_file)
                self.ligand.import_pdb("protein.pdb")
            except ImportError, e:
                sys.exit(1)

            if params.mode=="deform":
                self.structure_ligand=Protein()
                self.structure_ligand.import_pdb("protein.pdb")
                self.ligand.import_pdb("CA.pdb")

        else:
            #load monomeric structure
            self.ligand.import_pdb(params.ligand_file_name)


        #RECEPTOR STRUCTURE
        self.receptor = Protein()

        if params.receptor_style=="flexible":
            print ">> flexible docking requested for receptor, launching PCA..."
            try:
                self.flex_receptor=F.Flexibility_PCA()
                self.flex_receptor.compute_eigenvectors(params.receptor_topology,params.receptor_trajectory,params.receptor_align,params.receptor_ratio,params.mode,params.receptor_proj_file)
                self.receptor.import_pdb("protein.pdb")
            except ImportError, e:
                sys.exit(1)

            if params.mode=="deform":
                self.structure_receptor=Protein()
                self.structure_receptor.import_pdb("protein.pdb")
                self.receptor.import_pdb("CA.pdb")

        else:
            #load monomeric structure
            self.receptor.import_pdb(params.receptor_file_name)


        self.cg_atoms=[]

        if params.energy_type=="vdw":
            [self.index_ligand,self.index_receptor]=self.get_index(["CA"])
            #[self.index_ligand,self.index_receptor]=self.get_index(["CA","CB"])


    def get_index(self,atoms=["CA","CB"]):

        #generate a dummy assembly and extract the indexes where atoms of interest are located
        assembly = A.Assembly(self.ligand,self.receptor)
        assembly.place_ligand(np.array([0.0,0.0,0.0,0.0,0.0,0.0]))

        ligand_index=[]
        receptor_index=[]
        for aname in atoms:
        #append indexes of an element in atoms list for ligand
            [m,index]=assembly.atomselect_ligand("*","*",aname,True)
            for i in index:
                ligand_index.append(i)

            #append indexes of an element in atoms list for receptor
            [m,index]=assembly.atomselect_receptor("*","*",aname,True)
            for i in index:
                receptor_index.append(i)

        return [ligand_index,receptor_index]



class Space(S):
    def __init__(self,params,data):

        len_lig=0
        if params.ligand_style=="flexible":
            len_lig=len(data.flex_ligand.eigenspace_size)

        len_rec=0
        if params.receptor_style=="flexible":
            len_rec=len(data.flex_receptor.eigenspace_size)

        self.low=np.zeros(6+len_lig+len_rec)
        self.high=np.zeros(6+len_lig+len_rec)
        self.cell_size=np.zeros(6+len_lig+len_rec)
        self.boundary_type=np.zeros(6+len_lig+len_rec)


        #box size as given by receptor and ligand dimensions
        r_min=np.min(data.receptor.get_xyz(),axis=0)
        r_max=np.max(data.receptor.get_xyz(),axis=0)
        l_min=np.min(data.ligand.get_xyz(),axis=0)
        l_max=np.max(data.ligand.get_xyz(),axis=0)
        box_min=r_min-(l_max-l_min)
        box_max=r_max+(l_max-l_min)

        if len(params.high_input)!=len(params.low_input):
            print "ERROR: boundaryMin and boundaryMax should have the same length!"
            sys.exit(1)

        #assign low boundaries
        if params.low_input!="NA" :
            if len(params.low_input)==6:
                for i in xrange(0,len(params.low_input),1):
                    self.low[i]=params.low_input[i]
            else:
                print "ERROR: boundaryMin should contain 6 values (3 rotations, 3 translations)"
                sys.exit(1)
        else:
            print "WARNING: boundaryMin undefined, using default values"
            self.low[0]=box_min[0]
            self.low[1]=box_min[1]
            self.low[2]=box_min[2]
            self.low[3]=0.0
            self.low[4]=0.0
            self.low[5]=0.0

        #assign high boundaries
        if params.high_input!="NA" :
            if len(params.high_input)==6:
                for i in xrange(0,len(params.high_input),1):
                    self.high[i]=params.high_input[i]
            else:
                print "ERROR: boundaryMax should contain 6 values (3 rotation, 3 translation)"
                sys.exit(1)
        else:
            print "WARNING: boundaryMax undefined, using default values"
            self.high[0]=box_max[0]
            self.high[1]=box_max[1]
            self.high[2]=box_max[2]
            self.high[3]=360.0
            self.high[4]=360.0
            self.high[5]=360.0


        #add ligand eigenvector fluctuations in search space
        for i in xrange(0,len_lig,1):
            self.low[6+i]=-data.flex_ligand.eigenspace_size[i]
            self.high[6+i]=data.flex_ligand.eigenspace_size[i]

        #add receptor eigenvector fluctuations in search space
        for i in xrange(0,len_rec,1):
            self.low[6+len_lig+i]=-data.flex_receptor.eigenspace_size[i]
            self.high[6+len_lig+i]=data.flex_receptor.eigenspace_size[i]


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

        #check if target exists
        try: params.target
        except NameError:
            print 'ERROR: target measures not found'
            sys.exit(1)
        self.target=params.target
        
        self.constraint=params.constraint.split('.')[0]
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

        self.ligand_style=params.ligand_style
        self.receptor_style=params.receptor_style

        self.len_lig=0
        if params.ligand_style=="flexible":
            self.len_lig=len(self.data.flex_ligand.eigenspace_size)

        self.len_rec=0
        if params.receptor_style=="flexible":
            self.len_rec=len(self.data.flex_receptor.eigenspace_size)


    def evaluate(self,num,pos):

        exec 'import %s as constraint'%(self.constraint)
        import Assembly as A

        #if ligand is flexible, select the most appropriate frame
        if self.ligand_style=="flexible":
            deform_coeffs=pos[6:6+self.len_lig]
            if self.mode=="seed":
                pos_eig=self.data.flex_ligand.proj[:,self.data.flex_ligand.centroid]+deform_coeffs
                code,min_dist=vq(self.data.flex_ligand.proj.transpose(),np.array([pos_eig]))
                target_frame=min_dist.argmin()
                coords=self.data.flex_ligand.all_coords[:,target_frame]
                coords_reshaped=coords.reshape(len(coords)/3,3)
                self.data.ligand.set_xyz(coords_reshaped)
            else:
                coords=self.data.ligand.get_xyz()
                coords_reshaped=coords.reshape(len(coords)*3)

                for n in xrange(0,len(deform_coeffs),1):
                    coords_reshaped+=deform_coeffs[n]*self.data.flex_ligand.eigenvec[:,n]

                self.data.ligand.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))


        #if receptor is flexible, select the most appropriate frame
        if self.receptor_style=="flexible":
            deform_coeffs=pos[6+self.len_lig:6+self.len_lig++self.len_rec]
            if self.mode=="seed":
                pos_eig=self.data.flex_receptor.proj[:,self.data.flex_receptor.centroid]+deform_coeffs
                code,min_dist=vq(self.data.flex_receptor.proj.transpose(),np.array([pos_eig]))
                target_frame=min_dist.argmin()
                coords=self.data.flex_receptor.all_coords[:,target_frame]
                coords_reshaped=coords.reshape(len(coords)/3,3)
                self.data.receptor.set_xyz(coords_reshaped)
            else:
                coords=self.data.receptor.get_xyz()
                coords_reshaped=coords.reshape(len(coords)*3)

                for n in xrange(0,len(deform_coeffs),1):
                    coords_reshaped+=deform_coeffs[n]*self.data.flex_receptor.eigenvec[:,n]

                self.data.receptor.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))


        self.assembly = A.Assembly(self.data.ligand, self.data.receptor, self.data.cg_atoms)
        self.assembly.place_ligand(pos)

        #if needed, compute error with respect of target measures
        distance=0
        if len(self.target)!=0:

            measure = constraint.constraint_check(self.assembly)

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
        if len(self.data.index_ligand)>0:

            c1=0.2
            energy=self.interface_vdw()
            return c1*energy+(1-c1)*distance
            #return energy/len(self.data.index_ligand)+distance
        else:
            print "WHAT THE...???"
        #else:
        #   c1=0.001
        #   energy=self.measure_cg_energy(self.assembly,num)
        #   #fitness = coulomb+vdw+distance
        #   return c1*(energy[1]+energy[2])+(1-c1)*distance


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

        #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms
        m1=self.assembly.get_ligand_xyz()[self.data.index_ligand]
        m2=self.assembly.get_receptor_xyz()[self.data.index_receptor]

        #extract distances of every atom from all the others
        d=[]
        for i in xrange(0,len(m1),1):
            d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

        dist=np.array(d)

        #detect interfacing atoms (atom couples at less than a certain cutoff distance
        couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues

        energy=0
        for i in xrange(0,len(couples[0]),1):
            d=dist[couples[0,i],couples[1,i]]
            energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)

        return energy



class Postprocess(PP):

    def __init__(self,data,params):
        self.data=data
        self.params=params


        self.len_lig=0
        if params.ligand_style=="flexible":
            self.len_lig=len(self.data.flex_ligand.eigenspace_size)

        self.len_rec=0
        if params.receptor_style=="flexible":
            self.len_rec=len(self.data.flex_receptor.eigenspace_size)

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
    #threshold2 = clustering threshold
    def run(self) :

        exec 'import %s as constraint'%(self.constraint)

        #create output directory for generated PDB
        self.OUTPUT_DIRECTORY="result"
        if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
            os.mkdir(self.OUTPUT_DIRECTORY)

        clusters_file=open("%s/solutions.dat"%self.params.output_folder,"w")

        #use superclass method to filter acceptable solutions
        self.log=self.select_solutions(self.params)
        print ">> %s solutions filtered"%len(self.log)
        if len(self.log)==0:
            return

        #generate a dummy multimer and extract the indexes of C alphas for ligand and receptor
        multimer = A.Assembly(self.data.ligand, self.data.receptor, self.data.cg_atoms)
        multimer.place_ligand(np.array([0.0,0.0,0.0,0.0,0.0,0.0]))
        [m,index]=multimer.atomselect_ligand("*","*","CA",True)
        [m,index_ligand]=multimer.atomselect_ligand("*","*","CA",True)
        [m,index_receptor]=multimer.atomselect_receptor("*","*","CA",True)

        #load the monomeric structure positions
        s = Protein()
        if self.params.ligand_style=="flexible":
            s.import_pdb("protein.pdb")
        else:
            s.import_pdb(self.params.ligand_file_name)

        coords=s.get_xyz()

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

            #if ligand is flexible, select the most appropriate frame
            if self.params.ligand_style=="flexible":
                deform_coeffs=np.array(C[cnt-1])[6:6+self.len_lig]
                if self.params.mode=="seed":
                    pos_eig=self.data.flex_ligand.proj[:,self.data.flex_ligand.centroid]+deform_coeffs
                    code,min_dist=vq(self.data.flex_ligand.proj.transpose(),np.array([pos_eig]))
                    target_frame=min_dist.argmin()
                    coords=self.data.flex_ligand.all_coords[:,target_frame]
                    coords_reshaped=coords.reshape(len(coords)/3,3)
                    self.data.ligand.set_xyz(coords_reshaped)
                else:
                    coords=self.data.ligand.get_xyz()
                    coords_reshaped=coords.reshape(len(coords)*3)

                    for n in xrange(0,len(deform_coeffs),1):
                        coords_reshaped+=deform_coeffs[n]*self.data.flex_ligand.eigenvec[:,n]

                    self.data.ligand.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))


            #if receptor is flexible, select the most appropriate frame
            if self.params.receptor_style=="flexible":
                deform_coeffs=np.array(C[cnt-1])[6+self.len_lig:6+self.len_lig+self.len_rec]
                if self.params.mode=="seed":
                    pos_eig=self.data.flex_receptor.proj[:,self.data.flex_receptor.centroid]+deform_coeffs
                    code,min_dist=vq(self.data.flex_receptor.proj.transpose(),np.array([pos_eig]))
                    target_frame=min_dist.argmin()
                    coords=self.data.flex_receptor.all_coords[:,target_frame]
                    coords_reshaped=coords.reshape(len(coords)/3,3)
                    self.data.receptor.set_xyz(coords_reshaped)
                else:
                    coords=self.data.receptor.get_xyz()
                    coords_reshaped=coords.reshape(len(coords)*3)

                    for n in xrange(0,len(deform_coeffs),1):
                        coords_reshaped+=deform_coeffs[n]*self.data.flex_receptor.eigenvec[:,n]

                    self.data.receptor.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

            #create multimer
            pos = np.array(C[cnt-1])[0:6].astype(float)
            multimer1 = A.Assembly(self.data.ligand,self.data.receptor, self.data.cg_atoms)
            multimer1.place_ligand(pos)

            #write multimer
            multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))

            #clustering loop
            m1_1=multimer1.get_ligand_xyz()[index_ligand]
            m1_2=multimer1.get_receptor_xyz()[index_receptor]
            m1=np.concatenate((m1_1,m1_2),axis=0)

            cnt2=1
            for i in xrange(0,len(k),1) :

                #if ligand is flexible, select the most appropriate frame
                if self.params.ligand_style=="flexible":
                    deform_coeffs=np.array(P[k[i]])[6:6+self.len_lig]
                    if self.params.mode=="seed":
                        pos_eig=self.data.flex_ligand.proj[:,self.data.flex_ligand.centroid]+deform_coeffs
                        code,min_dist=vq(self.data.flex_ligand.proj.transpose(),np.array([pos_eig]))
                        target_frame1=min_dist.argmin()
                        coords=self.data.flex_ligand.all_coords[:,target_frame1]
                        coords_reshaped=coords.reshape(len(coords)/3,3)
                        self.data.ligand.set_xyz(coords_reshaped)
                    else:
                        coords=self.data.ligand.get_xyz()
                        coords_reshaped=coords.reshape(len(coords)*3)

                        for n in xrange(0,len(deform_coeffs),1):
                            coords_reshaped+=deform_coeffs[n]*self.data.flex_ligand.eigenvec[:,n]

                        self.data.ligand.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))


                #if receptor is flexible, select the most appropriate frame
                if self.params.receptor_style=="flexible":
                    deform_coeffs=np.array(P[k[i]])[6+self.len_lig:6+self.len_lig+self.len_rec]
                    if self.params.mode=="seed":
                        pos_eig=self.data.flex_receptor.proj[:,self.data.flex_receptor.centroid]+deform_coeffs
                        code,min_dist=vq(self.data.flex_receptor.proj.transpose(),np.array([pos_eig]))
                        target_frame2=min_dist.argmin()
                        coords=self.data.flex_receptor.all_coords[:,target_frame2]
                        coords_reshaped=coords.reshape(len(coords)/3,3)
                        self.data.receptor.set_xyz(coords_reshaped)
                    else:
                        coords=self.data.receptor.get_xyz()
                        coords_reshaped=coords.reshape(len(coords)*3)

                        for n in xrange(0,len(deform_coeffs),1):
                            coords_reshaped+=deform_coeffs[n]*self.data.flex_receptor.eigenvec[:,n]

                        self.data.receptor.set_xyz(coords_reshaped.reshape(len(coords_reshaped)/3,3))

                #extract positions of Calpha atoms in multimer
                multimer2 = A.Assembly(self.data.ligand,self.data.receptor,self.data.cg_atoms)
                multimer2.place_ligand(np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2],P[k[i]][3],P[k[i]][4],P[k[i]][5]]))
                m2_1=multimer1.get_ligand_xyz()[index_ligand]
                m2_2=multimer1.get_receptor_xyz()[index_receptor]
                m2=np.concatenate((m1_1,m1_2),axis=0)

                #compute RMSD within reference and current model
                rmsd=self.align(m1,m2)

                if rmsd<self.params.cluster_threshold :
                    cnt2+=1
                    P_C[k[i]]=cnt

            print ">>> clustered %s solutions on multimer %s"%(cnt2-1,cnt)

            #set centroid score with score of closest neighbor in set
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

        return

# Web site : http://lbm.epfl.ch


from Default import Parser as R
from Default import Space as S
from Default import Postprocess as PP

import numpy as np
import os, sys
from copy import deepcopy

from Protein import Protein
import Assembly as A


class Parser(R):

	def __init__(self):

		self.add('rmsdSelect','rmsd_select','array str',"NA")
                self.add('constraint','constraint','str',"NA")
                self.add('energy','energy_type','str',"vdw") 
                self.add('detectClash','detect_clash','str',"on")
		self.add('target','target','array float',np.array([]))
		self.add('ligand','ligand_file_name','str',"NA")
		self.add('receptor','receptor_file_name','str',"NA")

		self.add('cluster_threshold','cluster_threshold','float',"2")
		self.add('output_folder','output_folder','str',"result")

	def check_variables(self):

		if self.cluster_threshold<0:   
		    print "ERROR: clustering threshlod should be greater than 0!"
		    sys.exit(1)

		#check ligand files existence
		if self.ligand_file_name!="NA":
                    tmp=os.path.abspath(self.ligand_file_name)
		if os.path.isfile(self.ligand_file_name)!=1 :
        	    print "ERROR: ligand pdb file %s not found"%self.ligand_file_name
		    sys.exit(1)

		#check receptor existence
		if self.receptor_file_name!="NA":
                    tmp=os.path.abspath(self.receptor_file_name)
		if os.path.isfile(self.receptor_file_name)!=1 :
		    print "ERROR: receptor pdb file %s not found"%self.receptor_file_name
		    sys.exit(1)

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



class Data:
 
        index_ligand=[]
        index_receptor=[]
        cg_atoms=[]

	def __init__(self,params):
		self.ligand = Protein()
		self.ligand.import_pdb(params.ligand_file_name)

		self.receptor = Protein()
		self.receptor.import_pdb(params.receptor_file_name)

                self.cg_atoms=[]
 
		if params.energy_type=="vdw":
		    [self.index_ligand,self.index_receptor]=self.get_index(["CA","CB"])


	def get_index(self,atoms=["CA","CB"]):
	    
	    #generate a dummy assembly and extract the indexes where atoms of interest are located
	    assembly = A.Assembly(self.ligand,self.receptor)
	    assembly.place_ligand(np.array([0.0,0.0,0.0,0.0,0.0,0.0]))

	    ligand_index=[]
	    receptor_index=[]
	    for aname in atoms:
		#append indexes of an element in atoms list for ligand
		[m,index]=assembly.atomselect_ligand("*","*",aname,True)
		for i in index:
		    ligand_index.append(i)             

		#append indexes of an element in atoms list for receptor
		[m,index]=assembly.atomselect_receptor("*","*",aname,True)
		for i in index:
		    receptor_index.append(i)

	    return [ligand_index,receptor_index]



class Space(S):
	def __init__(self,params,data):


		self.low=np.zeros(6)
		self.high=np.zeros(6)
                #box size as given by receptor and ligand dimensions
		r_min=np.min(data.receptor.get_xyz(),axis=0)
		r_max=np.max(data.receptor.get_xyz(),axis=0)
		l_min=np.min(data.ligand.get_xyz(),axis=0)
		l_max=np.max(data.ligand.get_xyz(),axis=0)
		box_min=r_min-(l_max-l_min)
		box_max=r_max+(l_max-l_min)

		if len(params.high_input)!=len(params.low_input):
		    print "ERROR: boundaryMin and boundaryMax should have the same length!"
		    sys.exit(1)

		#assign low boundaries
		if params.low_input!="NA" :
		    if len(params.low_input)==6:
			for i in xrange(0,len(params.low_input),1):
			    self.low[i]=params.low_input[i]
		    else:
	     	        print "ERROR: boundaryMin should contain 6 values (3 rotations, 3 translations)"
			sys.exit(1)
		else:
		    print "WARNING: boundaryMin undefined, using default values"
		    self.low[0]=box_min[0]
		    self.low[1]=box_min[1]
		    self.low[2]=box_min[2]
		    self.low[3]=0.0
		    self.low[4]=0.0
		    self.low[5]=0.0

		#assign high boundaries
		if params.high_input!="NA" :
		    if len(params.high_input)==6:
			for i in xrange(0,len(params.high_input),1):
			    self.high[i]=params.high_input[i]
		    else:
		       print "ERROR: boundaryMax should contain 6 values (3 rotation, 3 translation)"
   		       sys.exit(1)
		else:
		    print "WARNING: boundaryMax undefined, using default values"
		    self.high[0]=box_max[0]
		    self.high[1]=box_max[1]
		    self.high[2]=box_max[2]
		    self.high[3]=360.0
		    self.high[4]=360.0
		    self.high[5]=360.0

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
		self.boundary_type=np.zeros(6)
		if params.boundary_type=="NA":
		    for i in xrange(0,3,1):
		        self.boundary_type[i]=0
		    for i in xrange(3,6,1):
		        self.boundary_type[i]=1
		elif params.boundary_type!="NA" and (np.logical_and(params.boundary_type!=0,params.boundary_type!=1)).any():
		    print 'ERROR: boundaries type not defined. This should be a list of 6 space-separated characters numbers (0=periodic, 1=repulsive)!'
		    sys.exit(1)
		else:
		    for i in xrange(0,6,1):
		        self.boundary_type[i]=params.boundary_type[i]



class Fitness:
	def __init__(self,data,params):

		#check if target exists
		try: params.target
		except NameError:
		    print 'ERROR: target measures not found'
		    sys.exit(1)
		self.target=params.target
		
		self.constraint=params.constraint.split('.')[0]
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
	

	def evaluate(self,num,pos):

            exec 'import %s as constraint'%(self.constraint)

            self.assembly = A.Assembly(self.data.ligand, self.data.receptor, self.data.cg_atoms)
            self.assembly.place_ligand(pos)

	    #if needed, compute error with respect of target measures
	    distance=0
	    if len(self.target)!=0:

	        measure = constraint.constraint_check(self.assembly)
	    	    
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
	    if len(self.data.index_ligand)>0:

                c1=0.2
                energy=self.interface_vdw()
                return c1*energy+(1-c1)*distance  
		#return energy/len(self.data.index_ligand)+distance
            else:
	        print "WHAT THE...???"
	    #else:
	    #	c1=0.001
	    #	energy=self.measure_cg_energy(self.assembly,num)  
	    #	#fitness = coulomb+vdw+distance  
	    #	return c1*(energy[1]+energy[2])+(1-c1)*distance


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

	    #extract coords of monomers 1 and 2 of multimeric structure according to desired atoms
	    m1=self.assembly.get_ligand_xyz()[self.data.index_ligand]
	    m2=self.assembly.get_receptor_xyz()[self.data.index_receptor]

	    #extract distances of every atom from all the others
	    d=[]
	    for i in xrange(0,len(m1),1):
	       d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))

	    dist=np.array(d)

	    #detect interfacing atoms (atom couples at less than a certain cutoff distance
	    couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
	    
	    energy=0
	    for i in xrange(0,len(couples[0]),1):    
		d=dist[couples[0,i],couples[1,i]]
		energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)

	    return energy



class Postprocess(PP):

	def __init__(self,data,params):
		self.data=data
		self.params=params

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
	#threshold2 = clustering threshold
	def run(self) :
	  
            exec 'import %s as constraint'%(self.constraint)

            #create output directory for generated PDB
            self.OUTPUT_DIRECTORY="result"
	    if os.path.isdir(self.OUTPUT_DIRECTORY)!=1:
		os.mkdir(self.OUTPUT_DIRECTORY)

	    clusters_file=open("%s/solutions.dat"%self.params.output_folder,"w")

	    #use superclass method to filter acceptable solutions
            self.log=self.select_solutions(self.params)
            print ">> %s solutions filtered"%len(self.log)
            if len(self.log)==0:
                return

	    #generate a dummy multimer and extract the indexes of C alpha
            multimer = A.Assembly(self.data.ligand, self.data.receptor, self.data.cg_atoms)
            multimer.place_ligand(np.array([0.0,0.0,0.0,0.0,0.0,0.0]))
            [m,index]=multimer.atomselect_ligand("*","*","CA",True)

            #load the monomeric structure positions
	    s = Protein()
	    s.import_pdb(self.params.ligand_file_name)   
            coords=s.get_xyz()

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

		#create multimer
		pos = np.array(C[cnt-1])[0:6].astype(float)
		multimer1 = A.Assembly(self.data.ligand,self.data.receptor, self.data.cg_atoms)
       		multimer1.place_ligand(pos)

		#write multimer
		multimer1.write_PDB("%s/assembly%s.pdb"%(self.OUTPUT_DIRECTORY,cnt))

		#clustering loop
                m1=multimer1.get_ligand_xyz()[index]
                cnt2=1
		for i in xrange(0,len(k),1) :
                    
                    self.data.ligand.set_xyz(coords)  
		    #multimer2 = A.Assembly(self.data.ligand)
		    multimer2 = A.Assembly(self.data.ligand,self.data.receptor,self.data.cg_atoms)
		    multimer2.place_ligand(np.array([P[k[i]][0],P[k[i]][1],P[k[i]][2],P[k[i]][3],P[k[i]][4],P[k[i]][5]]))
        	    m2=multimer2.get_ligand_xyz()[index]

		    rmsd=self.align(m1,m2)
                  		        
		    if rmsd<self.params.cluster_threshold :
                        cnt2+=1
		        P_C[k[i]]=cnt
		      
		print ">>> clustered %s solutions on multimer %s"%(cnt2,cnt)
  
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

		####generate output log###
		##write solution values
		#for item in C[cnt-1]:
		#    clusters_file.write("%s "%item)
		##write constraint values
		#for item in measure:
		#    clusters_file.write("%s "%item)
                ##write fitness
		#clusters_file.write("%s\n"%C_V[cnt-1])

	    return
