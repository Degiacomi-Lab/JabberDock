#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012-2014 EPFL (Ecole Polytechnique federale de Lausanne)
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
#hello
# Author : Matteo Degiacomi, matteo.degiacomi@gmail.com
# Web site : http://lbm.epfl.ch


import os, sys
import time
import numpy as np

from Default import Parser

from mpi4py import MPI
import dill
MPI.pickle.__init__(dill.dumps, dill.loads)
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

#check wether an appropriate number of parameters has been provided
if rank == 0:

    print('\n>         Parallel Optimization Workbench (POW)          <')
    print('> (c) 2012-14 Laboratory for Biomolecular Modeling, EPFL <')

    # checking if input is correct and has 3 arguments
    if len(sys.argv)!=2:
        print("\nERROR: expected 1 parameter as input!")
        print("USAGE: mpirun -n 4 POW.py input_file\n")
        sys.exit(1)

#get input file
infile=str(sys.argv[1])

if rank == 0 and os.path.isfile(infile)!=1 :
    print("ERROR: input file not found!")
    sys.exit(1)

#get program installation location and declare an environment variable (needed by some modules)
h=os.path.dirname(sys.argv[0])
home_dir=os.path.abspath(str(h))
os.environ['POW']=home_dir
sys.path.append(home_dir)

#also append subdirectories
sys.path.append("%s/modules"%home_dir)
sys.path.append("%s/optimizers"%home_dir)

#get local working directory and add to pythonpath
working_dir=os.getcwd()
sys.path.append(working_dir)

#parse input file, search for optimizer and module keywords
ptmp=Parser()
ptmp.add_standard()
ptmp.set_default_values()
ptmp.parse(infile,strict=False)

#add optimizer path in path, in case it's not in POW home
p,optm=os.path.split(ptmp.optimizer)
if len(p)>0:
    sys.path.append(p)

#add module path in path, in case it's not in POW home
p,mod=os.path.split(ptmp.module)
if len(p)>0:
    sys.path.append(p)

#preprocessing performed only by master node
#now necessary to be preprocessed by all due to pickling errors
#all prints by rank = 0
if rank > -1:

    if rank == 0:
        #init timer
        start=time.time()
        print('\n> PREPROCESSING\n')
        print(">> loading module %s..."%mod.split('.')[0])

    exec('import %s as mode'%(mod.split('.')[0]))
  
    #prepare input file parser 
    params=mode.Parser() 
    params.add_standard() 

    if rank == 0:
        print(">> loading optimizer %s..."%optm)

    exec('from %s import %s as Optimizer'%(optm.split('.')[0],optm.split('.')[0]))
    search=Optimizer()

    #if optimizer requires specific keywords, add them to parser
    if "add_keywords" in dir(search):
        params=search.add_keywords(params)

    #parser is ready, read input file and check values consistency
    if rank == 0:
        print('>> parsing input file...')
    params.set_default_values()

    params.parse(infile) 
    params.check_standard_variables() 
    params.check_variables() 

    #load requested data structures
    if rank == 0:
        print('>> importing data...')
    data=mode.Data(params)
   
    #build search space
    if rank == 0:
        print(">> generating search space...")

    space=mode.Space(params,data)

    if rank == 0:
        print(">> "+str(len(space.low))+"D space generated (min_pos, max_pos):")
        for i in range(0,len(space.low),1):
            print("   %s, %s"%(space.low[i],space.high[i]))

    #charge module fitness function or the one requested by user
    if params.fit=="NA":
        print(">> creating fitness function...")
        fitness=mode.Fitness(data,params)
    else:
        print(">>> importing user defined fitness function")
        try:
            exec('import %s as user_fit'%(params.fit.split('.')[0]))
        except ImportError as e:
            print("ERROR: load of user defined fitness function failed!")
            sys.exit(1)
        try:
            user_fit.Fitness
        except NameError:
            print('ERROR: user defined fitness function file should contain a class Fitness with a function evaluate(id,values_array) ')
        fitness=user_fit.Fitness(data,params)

else:
    params=None
    data=None
    space=None
    fitness=None
    search=None
    post=None

#propagate parameters, data, space, optimizer and fitness function to slaves
comm.Barrier()
#params=comm.bcast(params,root=0)
#space=comm.bcast(space,root=0)
#fitness=comm.bcast(fitness,root=0)
#data=comm.bcast(data,root=0)
#search=comm.bcast(search,root=0)
comm.Barrier()

#init optimization timer and launch optimization
if rank==0:
    start2=time.time()

search.launch(params,space,fitness)
comm.Barrier()

#launch postprocessing
if rank > -1:

    if rank == 0:
        end2=time.time()
        one_step=end2-start2
        print((">> optimization time (sec): {t:<8.3f}".format(t=one_step)))
        print('\n> POSTPROCESSING...')

    post=mode.Postprocess(data,params)
    
comm.Barrier()
#post=comm.bcast( post,root=0)
#comm.Barrier()
    
post.run()

if rank == 0:
    end=time.time()
    one_step=end-start
    print(("\n>> total execution time (sec): {t:<8.3f}\n".format(t=one_step)))
