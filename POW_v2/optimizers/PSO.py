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

import pickle
import shutil, os
import numpy as np
import copy
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

class PSO:

    def add_keywords(self,params):

        params.add('repeat','repeat','int',1)
        params.add('steps','max_steps','int',200)
        params.add('particles','n_particles','int',80)
        params.add('neighborSize','neigh_size','int',1)
        params.add('neighborType','neigh_type','str',"indexed")
        params.add('inertiaMax','inertia_max','float',0.9)
        params.add('inertiaMin','inertia_min','float',0.4)
        params.add('cp','cp','float',2.0)
        params.add('cn','cn','float',2.0)
        params.add('repulsion','repel','str',False)
        params.add('repulsion_factor','repel_factor', 'float', "NA")
        params.add('kar_threshold','kar','float',0.01)
        params.add('load_restart','restart_load','str',"NA")
        params.add('save_restart','restart_save','str',"swarm.restart")
        params.add('restart_freq','restart_freq','int',0)

        return params


    #just a function to check the optimizer-specific keywords passed via the input file
    def check_keywords(self):

        #if needed, init repellers
        if self.params.repel == 'on' :
            self.params.repel=True
            #repel_factor should be < 1!
            #self.repel_factor=0.05
        elif self.params.repel == 'off':
            self.params.repel=False
            self.params.repel_factor=0
        if self.params.repel != False and self.params.repel != True:
            print('ERROR: repulsion should either be "on" or "off"!')
            sys.exit(1)

        #check restart frequency
        if self.params.restart_freq>self.params.max_steps:
            print("ERROR: restart frequency should be smaller than total number of steps!")
            sys.exit(1)
        if self.params.restart_freq==0:
            self.params.restart_freq=int(float(self.params.max_steps)/10.0)
            #print ">>> a restart file will be created every %s steps"%self.restart_freq

        #check restart file existence
        if self.params.restart_load!="NA" and os.path.isfile(self.params.restart_load)!=1:
            print("ERROR: restart file %s not found"%self.params.restart_load)
            sys.exit(1)

        #check inertia values
        if self.params.inertia_max>1:
            print('ERROR: maximal inertia should be smaller than 1!')
            sys.exit(1)
        if self.params.inertia_min<0:
            print('ERROR: minimal inertia should be greater than 0!')
            sys.exit(1)
        if self.params.inertia_max<self.params.inertia_min :
            print('ERROR: maximal inertia is lower than minimal inertia (or one of these values are undefined)!')
            sys.exit(1)

        #check neighborhood conditions
        if self.params.neigh_type!="geographic" and self.params.neigh_type!="indexed":
            print("ERROR: neighborType should be either geographic either indexed")
            sys.exit(1)

        #check kick and reseed threshold value
        if self.params.kar<0:
            print("ERROR: kick and reseed (kar) threshold should be greater than 0!")
            sys.exit(1)


    def launch(self,params,space,fitness):

    #load parameterization, space, and fitness function information
        self.params=params
        self.space=space
        self.space.cell_size=self.space.high-self.space.low
       
        self.fitness=fitness

        if self.params.dimensions==-1:
            self.params.dimensions=len(self.space.low)

    #master node checks for optimizer-specific keywords consistency
        if rank ==0:
            self.check_keywords()

        #master node creates a swarm from a random state or restarted one
        if rank == 0:

            swarm=Swarm()

            #if restart file is provided, load the swarm state and prepare log file
            if self.params.restart_load!="NA":
                print(">> recovering swarm state from file %s..."%self.params.restart_load)
                if os.path.isfile(self.params.restart_load)==1 :
                    swarm.load_restart(self.params.restart_load)
                else:
                    print("ERROR: restart file %s not found!"%self.params.restart_load)
                    sys.exit(1)

                #if old logfile is available backup it and create a new one with the same name
                #containing the lines generated until the restart file time
                if os.path.isfile(self.params.output_file)==1:

                    shutil.copy2(self.params.output_file,"%s.bkp"%self.params.output_file)

                    fout_log = open(self.params.output_file, "w")
                    fin_log=open("%s.bkp"%self.params.output_file,"r")
                    totlines=self.params.n_particles*(swarm.repeat+1)*(swarm.ts+1)
                    for i in range(0,totlines,1):
                        line=fin_log.readline()
                        fout_log.write(line)

                    fin_log.close()
                    print(">> logging data in file %s, previous log backed up as %s.bkp"%(self.params.output_file, self.params.output_file))

                #if logfile is not found, create a new one
                else:
                    print("WARNING: previous log file %s not found, starting a new one...")
                    fout_log = open(self.params.output_file, "w")

                print(">> optimization restored from repetition %s and step %s"%(swarm.repeat+1,swarm.ts+1))

            #generate a random swarm otherwise
            else: #
                swarm.seed(self.params.n_particles,self.space) # indent in real, as part of else just above.
                fout_log = open(self.params.output_file, "w")  # indent in real


            fitness_best=1000000000

    # this is for the other processors
        else:
            swarm=None

        #synchronize all processes (get current timestep and repeat from swarm state)
        comm.Barrier()
        swarm=comm.bcast(swarm,root=0)
        start_repeat=swarm.repeat
        start_ts=swarm.ts
        comm.Barrier()


        #prepare repellers scaling factor, Xx the repeller should push according to the size of the cel right?
        if self.params.repel:
            self.params.scale_repeller=self.space.cell_size*self.params.repel_factor


        #################
        #PSO REPEAT LOOP#
        #################

        for r in range(start_repeat, self.params.repeat, 1):

            #init tracker for best fitness in the current repeat
            if rank == 0:
                print("\n> REPETITION %s\n"%(r+1))
                fitness_best=1000000000

            ###############
            #PSO MAIN LOOP#
            ###############

            for ts in range(start_ts,self.params.max_steps,1):

                #rescaling of inertia factor
                self.params.inertia=self.params.inertia_max-float(ts)/float(self.params.max_steps)*(self.params.inertia_max-self.params.inertia_min)

                #root spreads the data to everybody (itself included)
                comm.Barrier()
                swarm=comm.bcast(swarm,root=0)
                comm.Barrier()

                #every node elaborates the data. From received multidimensional array,
                #everybody extract the array corresponding to its rank,
                #and create a new one adding its rank to the extracted data
                p=rank
                list_pp=[]
                list_pv=[]
                list_pmp=[]
                list_pmv=[]
                list_cv=[]
                if self.params.repel:
                    list_repel=[]
                while p < self.params.n_particles :
                    # initialising particle
                    part=particle(p,self.fitness,self.space,self.params,swarm)
                    [part_pp,part_pv,part_pmp,part_pmv,part_cv, part_repel]=part.run()
                    #print "RANK %s of %s: %s -> %s"%(rank,size,p,part_cv)
                    list_pp.append(part_pp)
                    list_pv.append(part_pv)
                    list_pmp.append(part_pmp)
                    list_pmv.append(part_pmv)
                    list_cv.append(part_cv)
                    if self.params.repel and len(part_repel)>0:
                        list_repel.append(part_repel)
                    p+=size
                #print "RANK %s: data transmitted"%(rank)

                pp=np.array(list_pp)
                pv=np.array(list_pv)
                pmp=np.array(list_pmp)
                pmv=np.array(list_pmv)
                cv=np.array(list_cv)
                if self.params.repel:
                    rep=np.array(list_repel)

                #data updated by all CPU is sent back to root (root data included) in rank order,
                #in the form of a list
                comm.Barrier()
                pp_collect=comm.gather(pp,root=0)
                pv_collect=comm.gather(pv,root=0)
                pmp_collect=comm.gather(pmp,root=0)
                pmv_collect=comm.gather(pmv,root=0)
                cv_collect=comm.gather(cv,root=0)
                if self.params.repel:
                    r_collect=comm.gather(rep,root=0)

                #root will reshape the received lists in new ordered multidimensional arrays
                if rank == 0:
                    max_lines=int(np.ceil(float(self.params.n_particles)/float(size)))
                    swarm.part_pos=reformat(pp_collect,max_lines,size)
                    swarm.part_vel=reformat(pv_collect,max_lines,size)
                    swarm.part_max_pos=reformat(pmp_collect,max_lines,size)
                    swarm.part_max_val=reformat(pmv_collect,max_lines,size)
                    swarm.current_val=reformat(cv_collect,max_lines,size)

                    #concatenate new repelling points to the repellers list
                    if self.params.repel:
                        #iterate over processors summaries
                        cnt_rep=0
                        for j in range(0, len(r_collect), 1):
                            r_tmp=r_collect[j]
                            #iterate over list of repellers per processor
                            for i in range(0, len(r_tmp), 1):
                                cnt_rep+=1
                                #print ">>> now repelling from %s"%r_tmp[i]
                                swarm.repellers.append(r_tmp[i])

                        if cnt_rep>0:
                            print(">>> added %s new repellers (tot. repellers: %s)"%(cnt_rep,len(swarm.repellers)))

                    #write all positions to the log
                    for n in range(0,self.params.n_particles,1):
                        pos=swarm.part_pos[n]
                        #vel=swarm.part_vel[n]
                        fout_log.write("%s %s " % (r, n))
                        for item in pos:
                            fout_log.write("%s " % item)
                        fout_log.write("%s\n" % swarm.current_val[n])

                    if np.min(swarm.current_val) < fitness_best :
                        fitness_best=np.min(swarm.current_val)
                    print("step %s, best = %s"%(ts+1, fitness_best))

                    #if needed, write a swarm restart file (backup previous one, in case the computer crashes when writing the restart)
                    if self.params.restart_freq>0 and np.mod(ts+1, self.params.restart_freq)==0:
                        print(">> saving swarm state at repetition %s, step %s..."%(r+1,ts+1))
                        swarm.ts=ts
                        swarm.repeat=r
                        if os.path.isfile(self.params.restart_save)==1:
                            shutil.copy2(self.params.restart_save,"%s.old"%self.params.restart_save)
                        swarm.save_restart(self.params.restart_save)

            #finished a repetition, randomly reseed the swarm for the next one
            if rank == 0:
                swarm.seed(self.params.n_particles,self.space)



class Swarm:

    def __init__(self):
        #needed to build a proper restart (save ts and repeat at last restart time)
        self.ts=0
        self.repeat=0
        self.repellers=[]

    def seed(self,part,space):
        self.part_max_val = np.zeros(part)+10000
        self.current_val = np.zeros(part)+10000
        pp = []
        pv = []
        pmp = []
        for n in range(0,part,1):
            seed=rand_set(space.low,space.high)
            pp.append(seed)
            pmp.append(seed)
            pv.append(rand_set(-space.cell_size,space.cell_size))
        self.part_pos=np.array(pp)
        self.part_max_pos=np.array(pmp)
        self.part_vel=np.array(pv)


    #load swarm state from a restart file
    def load_restart(self,restart):
        file1=open(restart,"r")
        dataPickle=file1.read()
        file1.close()
        self.__dict__=pickle.loads(dataPickle)

    #save the swarm state into a restart file
    def save_restart(self,restart):
        file1=open(restart,"w")
        file1.write(pickle.dumps(self.__dict__))
        file1.close()



class particle :
    def __init__(self,n,fitness,space,params,swarm) :
        self.n=n
        self.fitness=fitness
        self.space=space
        self.params=params
        self.swarm=swarm

    def run(self) :

        #indicates wether to send a repeller position or not
        self.flag=False

        #copy locally the particle position and velocity
        self.row_pos = self.swarm.part_pos[self.n].copy()
        self.row_vel = self.swarm.part_vel[self.n].copy()

        #find best neighbor and save in self.best_pos the position of its best found fitness
        if self.params.neigh_type=="indexed":
            #indexed neighborhood
            if self.n==0:
                self.roi_val=np.array([self.swarm.part_max_val[self.params.n_particles-1],self.swarm.part_max_val[0],self.swarm.part_max_val[1]])
                self.roi_pos=np.array([self.swarm.part_max_pos[self.params.n_particles-1],self.swarm.part_max_pos[0],self.swarm.part_max_pos[1]])
            elif self.n==self.params.n_particles-1:
                self.roi_val=np.array([self.swarm.part_max_val[self.n-1],self.swarm.part_max_val[self.n],self.swarm.part_max_val[0]])
                self.roi_pos=np.array([self.swarm.part_max_pos[self.n-1],self.swarm.part_max_pos[self.n],self.swarm.part_max_pos[0]])
            else:
                self.roi_val=self.swarm.part_max_val[self.n-1:self.n+2]
                self.roi_pos=self.swarm.part_max_pos[self.n-1:self.n+2]
            self.best_value=np.min(self.roi_val)
            self.best_id=np.argmin(self.roi_val)
            self.best_pos=self.roi_pos[self.best_id]

        elif self.params.neigh_type=="geographic":
            #geographic neighborhood
            #particles distance from current particle (difference only, to gain time)
            self.dist=np.sum((self.swarm.part_pos-self.row_pos)**2,axis=1)
            #sort particles in order of distance
            #(current particle distance will always be zero, so we should exclude it)
            self.dist_sort=np.argsort(self.dist)
            #extract best fitness of neighbors (included particle's own fitness)
            self.best_id=self.swarm.part_max_val[self.dist_sort[0:self.params.neigh_size+1]].argmin()
            #extract best pos of particle having lowest best fitness within neighbors
            self.best_pos=self.swarm.part_pos[self.dist_sort[self.best_id]].copy()

        ###compute new velocity in search space###

        #inertia of previous speed
        self.previous=self.params.inertia*self.row_vel

        #influence of personal max
        self.local=self.swarm.part_max_pos[self.n]-self.row_pos
        #application perioding boundary conditions (for periodic dimensions)
        self.test=np.logical_and(abs(self.local)>abs(self.swarm.part_max_pos[self.n]+self.space.cell_size-self.row_pos),self.space.boundary_type==0)
        self.local[self.test]=self.swarm.part_max_pos[self.n][self.test]+self.space.cell_size[self.test]-self.row_pos[self.test]
        self.test=np.logical_and(abs(self.local)>abs(self.swarm.part_max_pos[self.n]-self.space.cell_size-self.row_pos),self.space.boundary_type==0)
        self.local[self.test]=self.swarm.part_max_pos[self.n][self.test]-self.space.cell_size[self.test]-self.row_pos[self.test]
        self.local_max=self.params.cn*np.random.rand(len(self.space.low))*(self.local)

        #influence of global max
        self.glob=self.best_pos-self.row_pos
        #application of periodic boundary conditions (for periodic dimensions)
        self.test=np.logical_and(abs(self.glob)>abs(self.best_pos+self.space.cell_size-self.row_pos),self.space.boundary_type==0)
        self.glob[self.test]=self.best_pos[self.test]+self.space.cell_size[self.test]-self.row_pos[self.test]
        self.test=np.logical_and(abs(self.glob)>abs(self.best_pos-self.space.cell_size-self.row_pos),self.space.boundary_type==0)
        self.glob[self.test]=self.best_pos[self.test]-self.space.cell_size[self.test]-self.row_pos[self.test]
        self.global_max=self.params.cp*np.random.rand(len(self.space.low))*(self.glob)

        #sum up contributions to get the new particle velocity
        self.row_vel_new_tmp=self.previous+self.local_max+self.global_max

        #rescale obtained velocity if particle is too fast
        self.test=np.abs(self.row_vel_new_tmp) > self.space.cell_size
        self.row_vel_new_tmp[self.test]=np.multiply(np.sign(self.row_vel_new_tmp[self.test]),self.space.cell_size[self.test])
        ### flag, kick, reseed ###
        #update particle position according to velocity if speed is high enough (wrapping on boundary conditions)
        self.velocity=float(np.sqrt(np.dot(self.row_vel_new_tmp,self.row_vel_new_tmp)))

        if self.velocity > self.params.kar :
            self.row_pos_new_tmp=self.row_pos+self.row_vel_new_tmp
            [self.row_pos_new,self.row_vel_new]=self.space.check_boundaries(self.row_pos_new_tmp,self.row_vel_new_tmp)
        #if velocity is too low, launch the kick and reseed procedure
        else :
            #if repellers are used, indicate that a flag should be placed
            if self.params.repel :
                self.flag=True

            #random kick particles being almost motionless and in a minima with bad fitness
            if self.swarm.current_val[self.n] > self.params.accept :
                #print ">>> kicking particle %s having fitness %s and velocity %s"%(self.n,self.swarm.current_val[self.n],self.velocity)
                print(">>> kicking particle %s (fitness = %s)"%(self.n,self.swarm.current_val[self.n]))
                self.row_vel_new_tmp=rand_set(-self.space.cell_size,self.space.cell_size)
                self.row_pos_new_tmp=self.row_pos+self.row_vel_new_tmp
                [self.row_pos_new,self.row_vel_new]=self.space.check_boundaries(self.row_pos_new_tmp,self.row_vel_new_tmp)

            #random kick and reseed particles being almost motionless and having a good fitness (needless to check)
            else :
                #print ">>> reseeding particle %s having fitness %s and velocity %s"%(self.n,self.swarm.current_val[self.n],self.velocity)
                print(">>> reseeding particle %s (fitness = %s)"%(self.n,self.swarm.current_val[self.n]))
                self.row_vel_new=rand_set(-self.space.cell_size,self.space.cell_size)
                self.row_pos_new=rand_set(self.space.low,self.space.high)
                self.swarm.part_max_val[self.n]=10000
                self.swarm.part_max_pos[self.n]=self.row_pos_new.copy()

        #influence of repulsion field
        if self.params.repel :
            self.bias=np.zeros(self.params.dimensions)
            for n in range(0,len(self.swarm.repellers),1):
                #if the particle is close to a repulsion point in a certain dimension, add its biasing contribution
                self.diff=self.row_pos-self.swarm.repellers[n]
                test = self.diff<self.space.cell_size*self.params.repel_factor
                if np.any(test):
                    #print "repelling now..."
                    self.bias[test]-=self.space.cell_size[test]*self.params.repel_factor/(self.diff[test]**2)
                    #print "done!"
            #bias velocity and position with repeller, and verify boundary conditions
            self.test=np.abs(self.bias) > self.space.cell_size
            self.bias_tmp=np.multiply(np.sign(self.bias[self.test]),self.space.cell_size[self.test])
            self.bias[self.test]=self.bias_tmp

            self.row_vel_new_tmp=self.row_vel_new+self.bias
            self.row_pos_new_tmp=self.row_pos_new+self.bias

            [self.row_pos_new,self.row_vel_new]=self.space.check_boundaries(self.row_pos_new_tmp,self.row_vel_new_tmp)
            #print "p%s: done particle %s"%(rank,self.n)

        self.f=self.fitness.evaluate(self.n,self.row_pos_new)

        #check whether a new max has been found and store all the max and max_pos in temporary buffers
        if self.f <= self.swarm.part_max_val[self.n] :
            self.max_pos_new=copy.copy(self.row_pos_new)
            self.max_val_new=copy.copy(self.f)
        else:
            self.max_pos_new=copy.copy(self.swarm.part_max_pos[self.n])
            self.max_val_new=copy.copy(self.swarm.part_max_val[self.n])

        #return  pos, vel, max pos, max fitness value and last fitness value
        if self.flag:
            return [self.row_pos_new, self.row_vel_new, self.max_pos_new, self.max_val_new, self.f, self.row_pos]
        else:
            return [self.row_pos_new, self.row_vel_new, self.max_pos_new, self.max_val_new, self.f, np.array([])]



def rand_set(min_a,max_a):
    out=np.zeros(len(min_a))
    for i in range(0,len(min_a),1):
        out[i]=min_a[i]+np.random.rand(1)[0]*(max_a[i]-min_a[i])
    return out

def reformat(product,max_lines,size):
    m=[]
    for i in range(0,max_lines,1):
        for s in range(0,size,1):
            tmp=product[s]
            try:
                l=tmp[i]
            except:
                l=None
            if l is not None:
                m.append(l)
    return np.array(m)
