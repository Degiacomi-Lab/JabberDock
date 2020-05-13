#!/usr/bin/env python

from __future__ import division  # future is >= 3.0, this code has mainly been used with 2.6 & 2.7
from __future__ import with_statement  # only necessary for python 2.5 and not in heavy use
# from __future__ import collections.MutableMapping # does not exist in future, otherwise 2.5 would work
from __future__ import print_function  # for cross-checking, available from python 2.6
import sys
if sys.version.startswith('3'):  # in python 3.x
    xrange = range
    raw_input = input

__version__ = "0.92.04 $Revision: 3322 $ $Date: 2012-11-22 18:05:10 +0100 (Thu, 22 Nov 2012) $"

import math
import random
import copy
import time  # not really essential
import collections, numpy as np # arange, cos, size, eye, inf, dot, floor, outer, zeros, linalg.eigh, sort, argsort, random, ones,...
from numpy import inf, array, dot, exp, log, sqrt, sum   # to access the built-in sum fct:  __builtins__.sum or del sum removes the imported sum and recovers the shadowed
#try:
#    import matplotlib.pylab as pylab  # also: use ipython -pylab
#    show = pylab.show
#    savefig = pylab.savefig   # we would like to be able to use cma.savefig() etc
#    closefig = pylab.close
#except:
#    pylab = None
#    print('  Could not import matplotlib.pylab, therefore ``cma.plot()`` etc. is not available')
#    def show():
#        pass

__docformat__ = "reStructuredText"  # this hides some comments entirely?

sys.py3kwarning = True  # TODO: out-comment from version 2.6



use_sent_solutions = True  # 5-30% CPU slower, particularly for large lambda, will be mandatory soon

# import other cma classes:
from FitnessFunctions import FitnessFunctions
from Misc import Misc
from NoiseHandler import NoiseHandler
from BestSolution import BestSolution
from BoundPenalty import BoundPenalty
from GenoPheno import GenoPheno
from CMAParameters import CMAParameters
from CMAViEParameters import CMAViEParameters
from CMAVIEConsParameters import CMAViEConsParameters
from DerivedDict import DerivedDictBase
from DerivedDict import SolutionDict
from ElapsedTime import ElapsedTime
from Rotation import Rotation
import os
# import bbobbenchmarks
import fgeneric

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

os.path.dirname(os.path.abspath(__file__))

class BlancClass(object):
    """blanc container class for having a collection of attributes"""

class OOOptimizer(object):

    def __init__(self, xstart, **more_args):
        """``xstart`` is a mandatory argument"""
        self.xstart = xstart
        self.more_args = more_args
        self.initialize()
    def initialize(self):
        """(re-)set to the initial state"""
        self.countiter = 0
        self.xcurrent = self.xstart[:]
        raise NotImplementedError('method initialize() must be implemented in derived class')
    def ask(self):
        """abstract method, AKA "get" or "sample_distribution", deliver new candidate solution(s), a list of "vectors"
        """
        raise NotImplementedError('method ask() must be implemented in derived class')
    def tell(self, solutions, function_values):
        """abstract method, AKA "update", prepare for next iteration"""
        self.countiter += 1
        raise NotImplementedError('method tell() must be implemented in derived class')
    def stop(self):
        """abstract method, return satisfied termination conditions in a dictionary like
        ``{'termination reason': value, ...}``, for example ``{'tolfun': 1e-12}``, or the empty
        dictionary ``{}``. The implementation of `stop()` should prevent an infinite loop.
        """
        raise NotImplementedError('method stop() is not implemented')
    def disp(self, modulo=None):
        """abstract method, display some iteration infos if ``self.iteration_counter % modulo == 0``"""
        raise NotImplementedError('method disp() is not implemented')
    def result(self):
        """abstract method, return ``(x, f(x), ...)``, that is, the minimizer, its function value, ..."""
        raise NotImplementedError('method result() is not implemented')

    def optimize(self, objectivefct, logger=None, verb_disp=100, iterations=None):

        if logger is None:
            if hasattr(self, 'logger'):
                logger = self.logger

        citer = 0
        while not self.stop():
            if iterations is not None and citer >= iterations:
                return self.result()
            citer += 1

            X = self.ask()         # deliver candidate solutions
            fitvals = [objectivefct(x) for x in X]

#            self.fileToWrite.write(str(citer)+"\t"+str( float(min(fitvals)) - 79.48)+"\n")


            self.tell(X, fitvals)  # all the work is done here

            self.disp(verb_disp)
            logger.add(self) if logger else None

        logger.add(self, modulo=bool(logger.modulo)) if logger else None
        if verb_disp:
            self.disp(1)
        if verb_disp in (1, True):
            print('termination by', self.stop())
            print('best f-value =', self.result()[1])
            print('solution =', self.result()[0])

        #self.fileToWrite.close()

        return self.result()[1]
    
    def optimize_CMA_POW(self, objectivefct, target, params=None, space=None, logger=None, verb_disp=100, iterations=None):

        if logger is None:
            if hasattr(self, 'logger'):
                logger = self.logger

        citer = 0
        fout_log = open(params.output_file, "a")
        logCMA = open("CMAlog.txt", "a")
                
        while not self.stop():
            if iterations is not None and citer >= iterations:
                return self.result()
            citer += 1
            
            # deliver candidate solutions
            X = self.ask_POW(space)         
#            X=self.ask()
            
            # extracting the PSO single fitness value (fitvas) and measures (3 rotations, 1 translation) as measures
            measures = []
            fitvals = []
            for x in X:
                fitvalstmp, measurestmp = objectivefct(0,x)
                measures.append( measurestmp )
                fitvals.append( fitvalstmp )

            
            # writing data to log file            
            for n in xrange(0,len(X),1):
                pos=copy.deepcopy(X[n])
                #vel=swarm.part_vel[n]
                fout_log.write("%s %s " % (0, n))
                logCMA.write("%s %s " % (0, n))
                
                # writing parameters
                for item in pos:
                    fout_log.write("%s " % item)
                    logCMA.write("%s " % item)
#                        fout_log.write("%s\n" % min(fitvals))
                
                fout_log.write("%s\n" % fitvals[n])
                
                logValues = [str(i) for i in measures[n]] # inversed the order of PSOfitness and logValues, bad code _> to be changed
                logCMA.write("%s " % (" ".join(logValues))) # 
                logCMA.write("%s\n" % fitvals[n])
            

            self.tell(X, fitvals)  # all the work is done here

            self.disp(verb_disp)
            logger.add(self) if logger else None

        logger.add(self, modulo=bool(logger.modulo)) if logger else None
        if verb_disp:
            self.disp(1)
        
        fout_log.close()
        logCMA.close()

        return self.result()[1]


    def optimize_ViE(self, objectivefct, target, maxfuneval, logger=None, verb_disp=100, iterations=None):
        '''This function was used for the BBOB benchmark, it is now use with POW protein assembly'''

        if logger is None:
            if hasattr(self, 'logger'):
                logger = self.logger

        citer = 0
        V = ViE()
        functions_eval = 0
        unsucess = 0
        maxFunctionEvals = 150000

        while not self.stop() and functions_eval < maxFunctionEvals:
            if iterations is not None and citer >= iterations:
                return self.result()
            citer += 1

            X = self.ask()         # deliver candidate solutions
#            fitvals = [objectivefct(0,x) for x in X]
            fitvals = [objectivefct(x) for x in X]

            

            functions_eval += len(X)

#            print ("===============================================================")
#            print (">> old mean: "+str(self.mean))
#            print ("step size: "+str(self.sigma))
#            print (">> evalauation of old mean: "+str(objectivefct(self.mean)))

            Xvie, fitvalsvie, survivorsNo = V.eleminateNonViable(X, fitvals)

            # update the CMA and ViE parameters only if there is at least one surviving solution
            if survivorsNo > 1:
#            if True:
                # updating CMA parameters, mu, sigma, ...
                self.tell_ViE(X, fitvals, survivorsNo, unsucess)  # all the work is done here

                # updating viability boundaries:
#                if survivorsNo == 0:
#                    unsucess += 1
#                elif len(Xvie) > 0:
                V.updateBoundaries(fitvalsvie)


            else:
                unsucess += 1

#                if unsucess % 100 == 0:
#                    print ("> Unsucessfuls iterations so far: "+str(unsucess))
#                    print ("fit: "+str(fitvals))
                if functions_eval >= maxFunctionEvals:
                    print ("___ function evals reached ___")


            self.disp(verb_disp)
            logger.add(self) if logger else None

        logger.add(self, modulo=bool(logger.modulo)) if logger else None
        if verb_disp:
            self.disp(1)
        if verb_disp in (1, True):
            print('termination by', self.stop())
            print('best f-value =', self.result()[1])
            print('solution =', self.result()[0])
#        print('best f-value =', self.result()[1])
#        #
#        print('\n __ best f-value =', self.result()[1])
#        print (">> final function eval: "+str(functions_eval))
        print('termination by', self.stop())
        print (">> unsuccessful generations: "+str(unsucess))
        print (">> unsuccessful functions eval: "+str(unsucess*12)+"\n")
        print ("")
        #print('solution =', self.result()[0])


        return self.result()[1]

    # ======================================================================================================================
    # ================================================ CMA VIE OPTIMIZATION ================================================
    # ======================================================================================================================

    def optimize_ViE_POW(self, objectivefct, target, params=None, space=None, logger=None, verb_disp=100, iterations=None):


        if rank == 0 and params != None:

            print (">> constraint optimization "+str(params.constrained_opt))

            citer = 0
            if params.constrained_opt:
                V = ViE_const(params)
                fout_log = open(params.output_file, "a")
                
            else:
                V = ViE()
                fout_log = open(params.output_file, "a")
            logCMA = open("CMAlog.txt", "a")


            unsucess = 0

            if logger is None:
                if hasattr(self, 'logger'):
                    logger = self.logger
                    
           


        citer = 0
        boundaryConvergedBool = 0 # this is used in the case of constrained optimization with AllAtom non bonded energy calculation
        maxFunctionEvals = params.evalByRestarts #8000 if params.constrained_opt else 3000 # was 4000 and 3000 
        functions_eval = 0
        termination = self.stop()

        constraintIterant = 0

        while not termination :

            if iterations is not None and citer >= iterations:
                return self.result()

            citer += 1

            # getting the population:
            if rank == 0:
                X = self.ask_POW(space,params)         # deliver candidate solutions
#                X = self.ask()


            else:
                X = None

            comm.Barrier()
            X = comm.bcast(X,root=0)
            comm.Barrier()


            # allocating jobs to cpus
            fitvalstmp=[]
            PSOfitnesstmp=[]
            cpu_no = rank
            cpu_order=[]

            while cpu_no < len(X):
#                print (X[cpu_no])
                if params.constrained_opt:
                    measures , fitness = objectivefct(boundaryConvergedBool,X[cpu_no])
                else:
                    measures, fitness = objectivefct(boundaryConvergedBool,X[cpu_no])                
                
                #print ("rank "+str(rank)+" energy: "+str(measures[-1]))
                
                fitvalstmp.append(measures)
                PSOfitnesstmp.append(fitness)

                cpu_order.append(cpu_no)
                cpu_no += size
#
            comm.Barrier()
            cpu_order_collect=comm.gather(cpu_order,root=0)
            fitvalstmp = comm.gather(fitvalstmp,root=0)
            PSOfitnesstmp= comm.gather(PSOfitnesstmp,root=0)
            comm.Barrier()

            if rank == 0:
                # extract right fitvals in the right order (the order in which the processors took the jobs) :

                ord1D = [v for x in xrange(0,len(cpu_order_collect),1) for v in cpu_order_collect[x] ]
                fitvals1D = [v for x in xrange(0,len(fitvalstmp),1) for v in fitvalstmp[x] ]
                fitvals = [0 for i in range(len(ord1D)) ]
                
                # also extract the single objective values of the PSO fitness function
                PSOfitness1D = [v for x in xrange(0,len(PSOfitnesstmp),1) for v in PSOfitnesstmp[x] ]
                PSOfitness = [0 for i in range(len(ord1D)) ]

                for i in xrange(0, len(ord1D), 1):
                    fitvals[ord1D[i]] = fitvals1D[i]
                    PSOfitness[ord1D[i]] = PSOfitness1D[i]
                    


                functions_eval += len(X)

                if 11 < 3:
                    print ("=================================")
                    print (X)
                    print(fitvals)
                    print ("=================================")

                # eliminate the individual solutions not fitting the boundary criteria
#                Xvie, fitvalsvie, survivorsNo = V.eleminateNonViable(X, fitvals)
                if params.constrained_opt:
                    if not params.updateSimultaneous:
                        killedIndex, survivorsIndex = V.eleminateNonViable(X, fitvals)
                    else:
                        killedIndex, survivorsIndex = V.eleminateNonViableSimultaneous(X, fitvals)

                    # formatting data to be passed to tell function
                    indeces = [killedIndex, survivorsIndex]
                    survivorsNo = len(survivorsIndex)
                    constraints= copy.deepcopy(fitvals)
                    fitvals = np.array(fitvals)[:,-1]

                    # extract the survivors coordinates a well as energy to be written in the log file in case of boundary convergence:
                    Xvie = [ X[i] for i in survivorsIndex  ]
                    fitvalsvie = [ fitvals[i] for i in survivorsIndex ]

                else:
                    Xvie, fitvalsvie, survivorsNo = V.eleminateNonViable(X, fitvals)

                    # if running POW, write the data in the log book:
                if not params.writeConverged :
                    
                    # first writing the protein parameters
                    for n in xrange(0,len(X),1):
                        pos=copy.deepcopy(X[n])
                        #vel=swarm.part_vel[n]
                        fout_log.write("%s %s " % (0, n))
                        logCMA.write("%s %s " % (0, n))
                        for item in pos:
                            fout_log.write("%s " % item)
                            logCMA.write("%s " % item)
#                        fout_log.write("%s\n" % min(fitvals))
                        
                        fout_log.write("%s\n" % fitvals[n])
                        
                        # then on a separate file, write the additional constraints and the PSO single fitness function
                        if params.constrained_opt:
                            logValues = [str(i) for i in constraints[n]]
                            logCMA.write("%s %s\n" % (" ".join(logValues),PSOfitness[n]))
                        else:
                            logValues = [str(i) for i in PSOfitness[n]] # inversed the order of PSOfitness and logValues, bad code _> to be changed
                            logCMA.write("%s " % (" ".join(logValues))) # 
                            logCMA.write("%s\n" % fitvals[n])


                # update the CMA and ViE parameters only if there is at least one surviving solution
                if survivorsNo > 0:

                    # updating CMA parameters, mu, sigma, ...
                    if params.constrained_opt:
                        self.tell_ViE_Constraint(X, fitvals, indeces, unsucess)

                        # change the X into the survivors to pass to the log file
                        X = Xvie
                        fitvals = fitvalsvie

                    else:
                        self.tell_ViE(X, fitvals, survivorsNo, unsucess)  # all the work is done here
                        V.updateBoundaries(fitvalsvie)

                    # GIORGIO -- increasing the sigma depending on the stopping criteria
#                    print (self.stopdict.xStop)
#                    if self.stopdict.xStop : #and not self.sigmaUpdated:
#                        print (">>> SIGMA UPDATE ON TERMINATION")
#                        self.sigma = 0.005
#                        self.sigmaUpdated = True


                else:
                    unsucess += 1

                    if unsucess % 100 == 0:
                        print ("> Unsucessfuls iterations so far: "+str(unsucess))
    #                    print ("fit: "+str(fitvals))
                    if functions_eval >= maxFunctionEvals:
                        print ("___ function evals reached ___")


                self.disp(verb_disp)
                logger.add(self) if logger else None

                # in this case stop only if the boundaries have converged
#                if not self.stop() and len(V.convergedBoundaryIndex) != (len(V.boundaries) -1) and functions_eval < maxFunctionEvals:

                # continue to update the energy while the boundaries are converged
                if not self.stop() and functions_eval < maxFunctionEvals:
                    termination = False
                else:
                    termination = True
                    print (">> Terminating CMA optimization process")



            # take the termination criteria from master node and propagate to slaves
            else:
                termination = None

            comm.Barrier()
            termination = comm.bcast(termination,root=0)
            comm.Barrier()

            if rank == 0 and params.constrained_opt and len(V.convergedBoundaryIndex) == (len(V.boundaries) -1):
                boundaryConvergedBool = 1
                
                  
                
            else:
                boundaryConvergedBool = 0
            comm.Barrier()
            boundaryConvergedBool = comm.bcast(boundaryConvergedBool,root=0)
            comm.Barrier()
                
            
            
            
        # -----  Termination (end of WHILE loop) -----

        if rank == 0:

            logger.add(self, modulo=bool(logger.modulo)) if logger else None
            if verb_disp:
                self.disp(1)
            if verb_disp in (1, True):
                print('termination by', self.stop())
                print('best f-value =', self.result()[1])


            print('termination by', self.stop())
            print (">> unsuccessful generations: "+str(unsucess))
            print (">> unsuccessful functions eval: "+str(unsucess*len(X))+"\n")
            print (">> total functions evaluations: "+str(functions_eval))
            print('solution =', self.result()[0])
            print ("")

            if params.constrained_opt:
                # only write to the file if the boundaries have converged(that you have a good solution)
                if params != None and len(V.convergedBoundaryIndex) == (len(V.boundaries) -1):
                    bool = 1
                    
                    
                    # in this case we're only writing the last best solutions which managed to solve all the constraints and minimize the energy
                    if params.writeConverged :
                        # first writing the protein parameters
                        for n in xrange(0,len(X),1):
                            pos=copy.deepcopy(X[n])
                            #vel=swarm.part_vel[n]
                            fout_log.write("%s %s " % (0, n))
                            logCMA.write("%s %s " % (0, n))
                            for item in pos:
                                fout_log.write("%s " % item)
                                logCMA.write("%s " % item)
    #                        fout_log.write("%s\n" % min(fitvals))
                            
                            fout_log.write("%s\n" % fitvals[n])
                            
                            # then write the additional constraints and the PSO single fitness function
                            if params.constrained_opt:
                                logValues = [str(i) for i in constraints[n]]
                                logCMA.write("%s %s\n" % (" ".join(logValues),PSOfitness[n]))
                            else:
                                logValues = [str(i) for i in PSOfitness[n]] # inversed the order of PSOfitness and logValues, bad code _> to be changed
                                logCMA.write("%s " % (" ".join(logValues))) # 
                                logCMA.write("%s\n" % fitvals[n])
                        

                else:
                    bool = 0
            else:
                bool = 0
            
            #closing the log files
            if params != None:
                fout_log.close()
                logCMA.close()


        # spreading the objects to return:
        else:
            bool = 0
            fitvals = []

        comm.Barrier()
        X = comm.bcast(X,root=0)
        fitvals = comm.bcast(fitvals, root=0)
        functions_eval = comm.bcast(functions_eval, root = 0)
        bool = comm.bcast(bool)
        comm.Barrier()

        return bool, functions_eval, X, fitvals

# ==================================================== VIABILITY EVOLUTION ==================================================== #




class ViE():
    def __init__(self, boundaries = None):
        self.boundaries = inf
        self.k=0.7
        self.updateBoundaryBool = True 

    def eleminateNonViable(self, X, fiteval):
        # extract the indexes of good solutions
        F = np.array(fiteval)
        totalInd = int(len(X))
        maxToKill = int(self.k * totalInd)
        print(F)
        # get only the X and F inside viability boundaries
        index = np.where(F <= self.boundaries )[0]
        Fviable = copy.deepcopy(F[index])
        Xviable = copy.deepcopy(np.array([X[i] for i in index]).reshape(len(index),len(X[0])))

        killedInd = totalInd - (len(index))
        # if the proportion to kill is not reached, kill until you do
        while killedInd < maxToKill:

            lowestFitIndex = np.where(Fviable == Fviable.max() )[0][0]
            Fviable = np.delete(Fviable, lowestFitIndex)
            Xviable = np.delete(Xviable, lowestFitIndex, axis=0)
            killedInd += 1

        survivors = totalInd - killedInd

        # check whether or not to update the boundary:
        if (killedInd <= maxToKill):
            self.updateBoundaryBool = True
        else:
            self.updateBoundaryBool = False



        return Xviable, Fviable, survivors

    def updateBoundaries(self, fiteval):
        if self.updateBoundaryBool:
            Fi = copy.deepcopy(np.array(fiteval))
            self.boundaries = Fi.max()
            print (">> new boundaries: "+str(self.boundaries)+"\n")


class ViE_const():
    def __init__(self, params):
        # store first the target boundaries -> boundary = [ [L(i),U(i)],...,[L(numberOfConstraints),U(numberOfConstraints)] ]
        self.targetBoundaries = [ [params.targetLow[x],params.targetHigh[x]] for x in xrange(0, len(params.targetLow),1) ]
        # the initial boundaries are relaxed
        self.boundaries = [[-inf,inf] for i in xrange(0,len(params.targetLow),1)] # -> [ [-inf,inf], [-inf, 18.3],..., -0.718 ]

        # add the boundary of the energy function:
        self.boundaries.append(inf) # __NEW__ _> was inf

        self.k=0.7
        self.updateBoundaryBool = True

        # this is used to iterate over the constraints, made it a class variable to remember the last index position of constraint array
        self.i = 0

        # this is used to know which boundaries are already converged:
        self.convergedBoundaryIndex = []

        # !!!! ONLY FOR CASE 3!!!!!, creation of a binary array to alternate between boundary updates:
        randBinList = lambda n: [random.randint(0,1) for b in range(1,n+1)]
        self.binaryArray = randBinList(len(self.boundaries)-1)


    def boundariesNotConverged(self):
        '''checks whether the target boundaries have been reached'''
        pass

    def eleminateNonViable(self, X, fitevals):
        '''just eliminates the individuals trespassing boundary values'''
        # define the maximal number of individuals to kill
        totalInd = int(len(X))
        maxToKill = int(self.k * totalInd)
        fiteval = np.array(fitevals)
        killedIndex = []

        #############################
        #print(self.boundaries)
        #print(np.shape(self.boundaries[1]))
        #print(fiteval)
        # Here we have a problem for JabberDock CMA POW. We only want one constraint, but we want to minimise another fitness function
        # Therefore we need to only loop through the dipole constraint (i.e. once), but CMA wants the Sc energy used as a constraint too
        # Therefore the following is hardcoded to only loop once - we can remove this if necessary by uncommenting
        ############################

        # instead of removing the individuals having constraints values beyond boundaries, their index is taken to remove them afterwards
        #for constno in xrange(0,len(self.boundaries),1):
        for constno in xrange(0,1,1):
            const = self.boundaries[constno]
            if type(const) is list:
                l = const[0]
                u = const[1]
            else: # for the energy
                l = -inf
                u = const

            for indNo in xrange (0,len(fiteval),1):
                #print('#########')
                #print(np.shape(fiteval), indNo, constno)
                #print(fiteval[indNo][constno])
                #print(fiteval[indNo][constno] < l) 
                #print(fiteval[indNo][constno] > u)
                #print(indNo not in killedIndex)
                if (fiteval[indNo][constno] < l or fiteval[indNo][constno] > u) and indNo not in killedIndex:
                    killedIndex.append(indNo)

        # extract indeces of survivors
        survivorsIndex = [v for v in xrange(0,len(fiteval),1) if v not in killedIndex]



        # in case the number of survivors is superior than the number of individuals to kill, eliminate until max number number to kill is reached
        while (len(killedIndex) < maxToKill):
            if self.i == (len(self.targetBoundaries) - len(self.convergedBoundaryIndex) ) or self.i == 0:
                self.i = 0
                self.randomIndeces = range(len(self.targetBoundaries))

                # remove the indeces which boundaries have already converged:
                for x in self.convergedBoundaryIndex:
                    self.randomIndeces.remove(x)

                random.shuffle(self.randomIndeces) # -> [0,2,3,1] for example


            if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1):
                c = self.boundaries[ self.randomIndeces[self.i] ] # -> [-inf, inf] i
                if type(c) is list: # in case of constraints
                    Ct = self.targetBoundaries[ self.randomIndeces[self.i] ] # -> [2, 3] i

            # sort the surviving individuals according to the current constraint (the one specified by the random index)
            ConstraintToIndexHash = {}
            nrjHash = {}
            for index in survivorsIndex:
                # extract the values of the current boundary that has to be updated
                if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1):
                    ConstraintToIndexHash[ fitevals[index][self.randomIndeces[self.i]] ] = index # -> { 18.888:0, 17.222:1, ... }
                # extract the energy values of the individuals
                nrjHash[ fitevals[index][-1] ] = index # -> { -12.23:0, -0.021:1, ... }

            # extract all the constraint values and sort them:
            constraintsOfSurvivors = copy.deepcopy( ConstraintToIndexHash.keys() )
            constraintsOfSurvivors.sort() # -> [1, 2, 3, 4, 18, 21, ...] i

            # extract the energy:
            energyOfSurvivors = copy.deepcopy( nrjHash.keys() )
            energyOfSurvivors.sort()


            # ---------------------------------------------------- CASE 1
            if 1 < 3:
                # case of the constraints
                print ("ENTERING UPDATING PART")
                if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1):
                    if ( (min(constraintsOfSurvivors) - c[0]) < (c[1] - max(constraintsOfSurvivors) ) or c[0] == -inf \
                         or ( max(constraintsOfSurvivors) < Ct[1] and min(constraintsOfSurvivors) < Ct[0] )  ) \
                         and min(constraintsOfSurvivors) <= Ct[0] and c[0] <= Ct[0] :
                        print("UPDATING LOWER BOUND")
                        # update the boundary
                        c[0] = min(constraintsOfSurvivors)

                        # send the index of the lowest individual so as to increase the number of killed individuals
                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[0]])

                        # kill the individual closest to the lower constraint/ remove the index from the suriving index
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[0]])
        #                constraintsOfSurvivors = constraintsOfSurvivors[1:]

                    # do the same as above but on the upper bounds
                    elif (c[1] >= Ct[1] and max(constraintsOfSurvivors) >= Ct[1]) :
                        print("UPDATING UPPER BOUND")
                        c[1] = constraintsOfSurvivors[-1]
                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
        #                constraintsOfSurvivors = constraintsOfSurvivors[:-1]

                    # checking if a boundary has converged
                    if min(constraintsOfSurvivors) >= Ct[0] and max(constraintsOfSurvivors) <= Ct[1]:
                        print ("BOUNDARY CONVERGED")
                        if self.randomIndeces[self.i] not in self.convergedBoundaryIndex:
                            self.convergedBoundaryIndex.append(self.randomIndeces[self.i])
                        print("state of converged index "+str(self.convergedBoundaryIndex))

                        # update the boundaries:
                        c[0] = Ct[0]
                        c[1] = Ct[1]
                        
                        self.i = -1

                    # with this you can see how each constraint as well their boundary values are updated
                    print ("constraints: "+str(constraintsOfSurvivors)+" "+str(self.randomIndeces[self.i]))
                    print (self.boundaries)
                    print ("----------------")

                # case of the energy
                else:
                    
                    
                    # this supposes that the energy function is always in the end of the boundaries array
                    self.boundaries[-1] = max(energyOfSurvivors)

                    # remove the individual having max energy
                    killedIndex.append(nrjHash[energyOfSurvivors[-1]])
                    survivorsIndex.remove(nrjHash[energyOfSurvivors[-1]])
                        
                    # with this you can see how each constraint as well their boundary values are updated
                    print ("Energy: "+str(energyOfSurvivors))
                    print (self.boundaries)
                    print ("----------------")



        # ---------------------------------------------------- CASE 2
            elif 11 < 3:
                # case of the constraints
                if type(c) is list:
                    if ( (min(constraintsOfSurvivors) - c[0]) < (c[1] - max(constraintsOfSurvivors) ) or c[0] == -inf) and min(constraintsOfSurvivors) < Ct[0] and c[0] < Ct[0] :


                        # send the index of the lowest individual so as to increase the number of killed individuals
                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[0]])

                        # kill the individual closest to the lower constraint/ remove the index from the suriving index
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[0]])
        #                constraintsOfSurvivors = constraintsOfSurvivors[1:]


                        # update the boundary not of the just killed individual but on the one after
                        delta = constraintsOfSurvivors[1] - c[0]
                        print ("delta "+str(delta))

                        c[0] = constraintsOfSurvivors[1] if constraintsOfSurvivors[1] < Ct[0] else constraintsOfSurvivors[0]

                        # update the upper boundary as well:
                        uc = max([Ct[1], (c[1] - delta)])
                        print ("u "+str(uc))

                        if uc != Ct[1]:
                            c[1] = uc

                            # eleminate the indeces higher than the new upper boundary:
                            for value in constraintsOfSurvivors:
                                if value > c[1]:
                                    killedIndex.append(ConstraintToIndexHash[value])
                                    survivorsIndex.remove(ConstraintToIndexHash[value])


                    # do the same as above but on the upper bounds
                    elif c[1] > Ct[1] and max(constraintsOfSurvivors) > Ct[1]:

                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])



                        # update the boundary not of the just killed individual but on the one after
                        print ("c[1] "+str(c[1])+" "+str(constraintsOfSurvivors[-2]))
                        delta =  c[1] - constraintsOfSurvivors[-2]
                        print ("delta "+str(delta))

                        c[1] = constraintsOfSurvivors[-2] if constraintsOfSurvivors[-2] > Ct[1] else constraintsOfSurvivors[-1]

                        # update the upper boundary as well:
                        lc = min([Ct[0], (c[0] - delta)])
                        print ("l "+str(lc))

                        if lc != Ct[0]:
                            c[0] = lc

                            # eleminate the indeces higher than the new upper boundary:
                            for value in constraintsOfSurvivors:
                                if value < c[0]:
                                    killedIndex.append(ConstraintToIndexHash[value])
                                    survivorsIndex.remove(ConstraintToIndexHash[value])

                    # checking if the boundary has converged
                    if min(constraintsOfSurvivors) > Ct[0] and max(constraintsOfSurvivors) < Ct[1]:
                        print ("BOUNDARY CONVERGED")
                        if self.randomIndeces[self.i] not in self.convergedBoundaryIndex:
                            self.convergedBoundaryIndex.append(self.randomIndeces[self.i])
                        print("state of converged index "+str(self.convergedBoundaryIndex))

                        # update the boundaries:
                        c[0] = Ct[0]
                        c[1] = Ct[1]

                        self.i = -1

                # case of the energy
                else:
                    # this supposes that the energy function is always in the end of the boundaries array
                    self.boundaries[-1] = max(constraintsOfSurvivors)

                    killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                    survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])

                # with this you can see how each constraint as well their boundary values are updated
                print ("constraints: "+str(constraintsOfSurvivors)+" "+str(self.randomIndeces[self.i]))
                print (self.boundaries)
                print ("----------------")

            # ---------------------------------------------------- CASE 3
            elif 11 < 3:
                # case of the constraints
                if type(c) is list:
                    if min(constraintsOfSurvivors) < Ct[0] and c[0] < Ct[0] and self.binaryArray[self.randomIndeces[self.i]] == 0 :
                        print ("Binary array before: "+str(self.binaryArray))
                        # update the boundary
                        c[0] = min(constraintsOfSurvivors)

                        # send the index of the lowest individual so as to increase the number of killed individuals
                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[0]])

                        # kill the individual closest to the lower constraint/ remove the index from the suriving index
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[0]])
        #                constraintsOfSurvivors = constraintsOfSurvivors[1:]

                        # change the number of the binary array in order to alternate to the next boundary:
                        self.binaryArray[self.randomIndeces[self.i]] = 1
                        print ("Binary array NOW: "+str(self.binaryArray))

                    # do the same as above but on the upper bounds
                    elif c[1] > Ct[1] and max(constraintsOfSurvivors) > Ct[1] and self.binaryArray[self.randomIndeces[self.i]] == 1:
                        print ("Binary array before: "+str(self.binaryArray))
                        c[1] = constraintsOfSurvivors[-1]
                        killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                        survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                        self.binaryArray[self.randomIndeces[self.i]] = 0
                        print ("Binary array NOW: "+str(self.binaryArray))

                    # checking if the boundary has converged
                    if min(constraintsOfSurvivors) > Ct[0] and max(constraintsOfSurvivors) < Ct[1]:
                        print ("BOUNDARY CONVERGED")
                        if self.randomIndeces[self.i] not in self.convergedBoundaryIndex:
                            self.convergedBoundaryIndex.append(self.randomIndeces[self.i])
                        print("state of converged index "+str(self.convergedBoundaryIndex))

                        # update the boundaries:
                        c[0] = Ct[0]
                        c[1] = Ct[1]

                        self.i = -1

                # case of the energy
                else:
                    # this supposes that the energy function is always in the end of the boundaries array
                    self.boundaries[-1] = max(constraintsOfSurvivors)

                    killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                    survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])



                # with this you can see how each constraint as well their boundary values are updated
                print ("constraints: "+str(constraintsOfSurvivors)+" "+str(self.randomIndeces[self.i]))
                print (self.boundaries)
                print ("----------------")

            self.i += 1

#        print ("EXITED WHILE LOOP _-----------")
#        print (killedIndex)
#        print(">> Survivors: "+str(survivorsIndex))

        return killedIndex, survivorsIndex













# ==================================================== COVARIANCE MATRIX ADAPTATION EVOLUTION STRATEGY ====================================================

    def eleminateNonViableSimultaneous(self, X, fitevals):
        '''just eliminates the individuals trespassing boundary values'''
        # define the maximal number of individuals to kill
        totalInd = int(len(X))
        maxToKill = int(self.k * totalInd)
        fiteval = np.array(fitevals)
        killedIndex = []


        # instead of removing the individuals having constraints values beyond boundaries, their index is taken to remove them afterwards
        for constno in xrange(0,len(self.boundaries),1):
            const = self.boundaries[constno]
            if type(const) is list:
                l = const[0]
                u = const[1]
            else: # for the energy
                l = -inf
                u = const

            for indNo in xrange (0,len(fiteval),1):
                if (fiteval[indNo][constno] < l or fiteval[indNo][constno] > u) and indNo not in killedIndex:
                    killedIndex.append(indNo)

        # extract indeces of survivors
        survivorsIndex = [v for v in xrange(0,len(fiteval),1) if v not in killedIndex]



        # in case the number of survivors is superior than the number of individuals to kill, eliminate until max number to kill is reached
        while (len(killedIndex) < maxToKill):
            # conditions to reset i, the iterator of random indeces
#             if self.i == (len(self.targetBoundaries) - len(self.convergedBoundaryIndex) ) or self.i == 0:
            if self.i == 0 or self.i == (len(self.boundaries) - len(self.convergedBoundaryIndex)):
                
                self.i = 0
                if 11 < 3:
                    self.randomIndeces = range(len(self.targetBoundaries))
                else:
                    self.randomIndeces = range(len(self.boundaries) )

                # remove the indeces which boundaries have already converged:
                for x in self.convergedBoundaryIndex:
                    self.randomIndeces.remove(x)

                random.shuffle(self.randomIndeces) # -> [0,2,3,1] for example

            if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1):
                c = self.boundaries[ self.randomIndeces[self.i] ] # -> [-inf, inf] i
                if type(c) is list: # in case of geometric constraints
                    Ct = self.targetBoundaries[ self.randomIndeces[self.i] ] # -> [2, 3] i
            
            # ---------------------- GETTING INDECES AND VALUES OF SURVIVORS
            
            # sort the surviving individuals according to the current constraint (the one specified by the random index)
            ConstraintToIndexHash = {}
            nrjHash = {}
            for index in survivorsIndex:
                # extract the values of the current boundary that has to be updated
                if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1):
                    ConstraintToIndexHash[ fitevals[index][self.randomIndeces[self.i]] ] = index # -> { 18.888:0, 17.222:1, ... }
                # extract the energy values of the individuals
                nrjHash[ fitevals[index][-1] ] = index # -> { -12.23:0, -0.021:1, ... }

            # extract all the constraint values and sort them:
            constraintsOfSurvivors = copy.deepcopy( ConstraintToIndexHash.keys() )
            constraintsOfSurvivors.sort() # -> [1, 2, 3, 4, 18, 21, ...] i

            # extract the energy:
            energyOfSurvivors = copy.deepcopy( nrjHash.keys() )
            energyOfSurvivors.sort()


            # ---------------------------------------------------- CASE 1
#             if 1 < 3:
            # case of the constraints
            print ("> ENTERING UPDATING PART")
            if len(self.convergedBoundaryIndex) != (len(self.boundaries) -1) and type(c) is list:
#                 if type(c) is list:
                if ( (min(constraintsOfSurvivors) - c[0]) < (c[1] - max(constraintsOfSurvivors) ) or c[0] == -inf \
                     or ( max(constraintsOfSurvivors) < Ct[1] and min(constraintsOfSurvivors) < Ct[0] )  ) \
                     and min(constraintsOfSurvivors) <= Ct[0] and c[0] <= Ct[0] :
                    print("> UPDATING LOWER BOUND")
                    # update the boundary
                    c[0] = min(constraintsOfSurvivors)

                    # send the index of the lowest individual so as to increase the number of killed individuals
                    killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[0]])

                    # kill the individual closest to the lower constraint/ remove the index from the suriving index
                    survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[0]])
    #                constraintsOfSurvivors = constraintsOfSurvivors[1:]

                # do the same as above but on the upper bounds
                elif (c[1] >= Ct[1] and max(constraintsOfSurvivors) >= Ct[1]) :
                    print("> UPDATING UPPER BOUND")
                    c[1] = constraintsOfSurvivors[-1]
                    killedIndex.append(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
                    survivorsIndex.remove(ConstraintToIndexHash[constraintsOfSurvivors[-1]])
    #                constraintsOfSurvivors = constraintsOfSurvivors[:-1]

                # checking if a boundary has converged
                if min(constraintsOfSurvivors) >= Ct[0] and max(constraintsOfSurvivors) <= Ct[1]:
                    print ("> BOUNDARY CONVERGED <")
                    if self.randomIndeces[self.i] not in self.convergedBoundaryIndex:
                        self.convergedBoundaryIndex.append(self.randomIndeces[self.i])
                    print("state of converged index "+str(self.convergedBoundaryIndex))

                    # update the boundaries:
                    c[0] = Ct[0]
                    c[1] = Ct[1]
                    
                    self.i = -1

                # with this you can see how each constraint as well their boundary values are updated
                print ("constraints: "+str(constraintsOfSurvivors)+" "+str(self.randomIndeces[self.i]))
                print (self.boundaries)
                print ("----------------")

            # case of the energy
            else:
                
                # this supposes that the energy function is always in the end of the boundaries array
                self.boundaries[-1] = max(energyOfSurvivors)

                # remove the individual having max energy
                killedIndex.append(nrjHash[energyOfSurvivors[-1]])
                survivorsIndex.remove(nrjHash[energyOfSurvivors[-1]])
                    
                # with this you can see how each constraint as well their boundary values are updated
                print ("> Energy: "+str(energyOfSurvivors))
                print (self.boundaries)
                print ("----------------")

            self.i += 1

#        print ("EXITED WHILE LOOP _-----------")
#        print (killedIndex)
#        print(">> Survivors: "+str(survivorsIndex))

        return killedIndex, survivorsIndex


class CMAEvolutionStrategy(OOOptimizer):

    # __all__ = ()  # TODO this would be the interface

    #____________________________________________________________
    @property  # read only attribute decorator for a method
    def popsize(self):
        """number of samples by default returned by` ask()`
        """
        return self.sp.popsize

    # this is not compatible with python2.5:
    #     @popsize.setter
    #     def popsize(self, p):
    #         """popsize cannot be set (this might change in future)
    #         """
    #         raise _Error("popsize cannot be changed (this might change in future)")

    #____________________________________________________________
    #____________________________________________________________
    def stop(self, check=True):
        """return a dictionary with the termination status.
        With ``check==False``, the termination conditions are not checked and
        the status might not reflect the current situation.
        """

        if (check and self.countiter > 0 and self.opts['termination_callback'] and
                self.opts['termination_callback'] != str(self.opts['termination_callback'])):
            self.callbackstop = self.opts['termination_callback'](self)

        return self.stopdict(self if check else None)  # update the stopdict and return a Dict

    #____________________________________________________________
    #____________________________________________________________
    def __init__(self, x0, sigma0, inopts = {}):
        """see class `CMAEvolutionStrategy`

        """
#        # writing results to a file:
#        self.fileToWrite = open("f1_CMAVIE.txt","w")

        self.inputargs = dict(locals()) # for the record
        del self.inputargs['self'] # otherwise the instance self has a cyclic reference
        self.inopts = inopts
        opts = Options(inopts).complement()  # Options() == fmin([],[]) == defaultOptions()

        if opts['noise_handling'] and eval(opts['noise_handling']):
            raise ValueError('noise_handling not available with class CMAEvolutionStrategy, use function fmin')
        if opts['restarts'] and eval(opts['restarts']):
            raise ValueError('restarts not available with class CMAEvolutionStrategy, use function fmin')


        # ---------------------------- DEFINING THE MEAN ---------------------------
        if x0 == str(x0):
            x0 = eval(x0)
        self.mean = array(x0)  # should not have column or row, is just 1-D
#        print ("self.mean:")
#        print (self.mean)
        if self.mean.ndim == 2:
            print('WARNING: input x0 should be a list or 1-D array, trying to flatten ' +
                    str(self.mean.shape) + '-array')
            if self.mean.shape[0] == 1:
                self.mean = self.mean[0]
            elif self.mean.shape[1] == 1:
                self.mean = array([x[0] for x in self.mean])
        if self.mean.ndim != 1:
            raise _Error('x0 must be 1-D array')
        if len(self.mean) <= 1:
            raise _Error('optimization in 1-D is not supported (code was never tested)')

        self.N = self.mean.shape[0]
        N = self.N
        self.mean.resize(N) # 1-D array, not really necessary?!
        self.x0 = self.mean
        self.mean = self.x0.copy()  # goes to initialize

        # ---------------------------- DEFINING SIGMA
        self.sigma0 = sigma0
        if isinstance(sigma0, str):  # TODO: no real need here (do rather in fmin)
            self.sigma0 = eval(sigma0)  # like '1./N' or 'np.random.rand(1)[0]+1e-2'
        if np.size(self.sigma0) != 1 or np.shape(self.sigma0):
            raise _Error('input argument sigma0 must be (or evaluate to) a scalar')
        self.sigma = self.sigma0  # goes to inialize

        # ---------------------------- DEFINING/EXPANDING OPTIONS
        # extract/expand options
        opts.evalall(locals())  # using only N
        self.opts = opts

        self.randn = opts['randn']
        self.gp = GenoPheno(N, opts['scaling_of_variables'], opts['typical_x'],
            opts['bounds'], opts['fixed_variables'], opts['transformation'])
        self.boundPenalty = BoundPenalty(self.gp.bounds)
        s = self.gp.geno(self.mean)
        self.mean = self.gp.geno(self.mean, bounds=self.gp.bounds)
        self.N = len(self.mean)
#        print ("self.N:")
#        print (self.N)
        N = self.N
        if (self.mean != s).any():
            print('WARNING: initial solution is out of the domain boundaries:')
            print('  x0   = ' + str(self.inputargs['x0']))
            print('  ldom = ' + str(self.gp.bounds[0]))
            print('  udom = ' + str(self.gp.bounds[1]))
        self.fmean = np.NaN             # TODO name should change? prints nan (OK with matlab&octave)
        self.fmean_noise_free = 0.  # for output only

        # ---------------------------- DEFINING POPULATION
        self.sp = CMAParameters(N, opts)
#        print ("self.sp.popsize: "+str(self.sp.popsize))
        self.sp0 = self.sp  # looks useless, as it is not a copy

        # initialization of state variables
        self.countiter = 0
        self.countevals = max((0, opts['verb_append'])) if type(opts['verb_append']) is not bool else 0
        self.ps = np.zeros(N)
        self.pc = np.zeros(N)

        stds = np.ones(N)
        self.sigma_vec = np.ones(N) if np.isfinite(self.sp.dampsvec) else 1
        if np.all(self.opts['CMA_teststds']):  # also 0 would not make sense
            stds = self.opts['CMA_teststds']
            if np.size(stds) != N:
                raise _Error('CMA_teststds option must have dimension = ' + str(N))
        if self.opts['CMA_diagonal']:  # is True or > 0
            # linear time and space complexity
            self.B = array(1) # works fine with np.dot(self.B, anything) and self.B.T
            self.C = stds**2  # TODO: remove this!?
            self.dC = self.C
        else:
            self.B = np.eye(N) # identity(N), do not from matlib import *, as eye is a matrix there
            # prevent equal eigenvals, a hack for np.linalg:
            self.C = np.diag(stds**2 * exp(1e-6*(np.random.rand(N)-0.5)))
            self.dC = np.diag(self.C)
            self.Zneg = np.zeros((N, N))
        self.D = stds

        self.flgtelldone = True
        self.itereigenupdated = self.countiter
        self.noiseS = 0  # noise "signal"
        self.hsiglist = []

        if not opts['seed']:
            np.random.seed()
            six_decimals = (time.time() - 1e6 * (time.time() // 1e6))
            opts['seed'] = 1e5 * np.random.rand() + six_decimals + 1e5 * (time.time() % 1)
        opts['seed'] = int(opts['seed'])
        np.random.seed(opts['seed'])

        # ---------------------------- CREATING SOLUTION DICTIONNARIES
        self.sent_solutions = SolutionDict()
        self.best = BestSolution()

        out = {}  # TODO: obsolete, replaced by method results()?
        out['best'] = self.best
        # out['hsigcount'] = 0
        out['termination'] = {}
        self.out = out

        self.const = BlancClass()
        self.const.chiN = N**0.5*(1-1./(4.*N)+1./(21.*N**2)) # expectation of norm(randn(N,1))

        # ---------------------------- DEFINING STOPPING CRITERIA
        # attribute for stopping criteria in function stop
        self.stopdict = CMAStopDict()
        self.callbackstop = 0

        self.fit = BlancClass()
        self.fit.fit = []   # not really necessary
        self.fit.hist = []  # short history of best
        self.fit.histbest = []   # long history of best
        self.fit.histmedian = [] # long history of median

        self.more_to_write = []  #[1, 1, 1, 1]  #  N*[1]  # needed when writing takes place before setting

        # say hello
        if opts['verb_disp'] > 0:
            sweighted = '_w' if self.sp.mu > 1 else ''
            smirr = 'mirr%d' % (self.sp.lam_mirr) if self.sp.lam_mirr else ''
#            print('(%d' % (self.sp.mu) + sweighted + ',%d' % (self.sp.popsize) + smirr + ')-CMA-ES' +
#                  ' (mu_w=%2.1f,w_1=%d%%)' % (self.sp.mueff, int(100*self.sp.weights[0])) +
#                  ' in dimension %d (seed=%d, %s)' % (N, opts['seed'], time.asctime())) # + func.__name__
            if opts['CMA_diagonal'] and self.sp.CMA_on:
                s = ''
                if opts['CMA_diagonal'] is not True:
                    s = ' for '
                    if opts['CMA_diagonal'] < np.inf:
                        s += str(int(opts['CMA_diagonal']))
                    else:
                        s += str(np.floor(opts['CMA_diagonal']))
                    s += ' iterations'
                    s += ' (1/ccov=' + str(round(1./(self.sp.c1+self.sp.cmu))) + ')'
                print('   Covariance matrix is diagonal' + s)

    # ---------------------------- FUNCTION TO ASK GENO PHENO POW MODE! --------------------

    def ask_POW(self, space, params, number=None, xmean=None, sigma_fac=1):

        if 1 < 3:
            pop_geno = self.ask_geno(number, xmean, sigma_fac) 
        else:
            pop_geno = self.ask_geno_resampled(space, number, xmean, sigma_fac)

        # variable used to print ot not the current phenotype and genotype
        printVar = 0

        if printVar:
            print ("GENOTYPE")
            for i in pop_geno:
                print (i)
            

        # transform the genotype depending on the boundary conditions:
        v = copy.deepcopy(pop_geno)
        
        np.savetxt("CMAVIE_coords.dat", v)
#         strLast=[str(l) for l in params.last]
#         params.CMAVIE_coords.write(strLast+"\n")

        # transforming the solutions parameters in the orginal range of 0 to 1 to 0-360 * 3 and then translation value range (eg. 2-6)
        pop_phenoTemp = [ ((x*(space.high-space.low))+space.low)  for x in pop_geno]
        

        if printVar:
            print ("PHENO TEMP")
            for i in pop_phenoTemp:
                print (i)

        pop_pheno = []
        
        # checking for the periodic boundaries consistency (POW)
        for x in pop_phenoTemp:
            xi, vi = space.check_boundaries(x,v[0])
            pop_pheno.append(xi)

        if printVar:
            print ("pop_pheno")
            for i in pop_pheno:
                print (i)


        if len(pop_pheno) != len(pop_geno):
            print ("pheno-geno not equal")

        if not self.gp.isidentity or use_sent_solutions:  # costs 25% in CPU performance with N,lambda=20,200
            # archive returned solutions, first clean up archive
            if self.countiter % 30/self.popsize**0.5 < 1:
                self.sent_solutions.truncate(0, self.countiter - 1 - 3 * self.N/self.popsize**0.5)
            # insert solutions
            for i in xrange(len(pop_geno)):
                self.sent_solutions[pop_pheno[i]] = {'geno': pop_geno[i],
                                            'pheno': pop_pheno[i],
                                            'iteration': self.countiter}



        return pop_pheno



    # ---------------------------- FUNCTION TO ASK GENO -> PHENO ---------------------------
    def ask_geno_resampled(self, space, number=None, xmean=None, sigma_fac=1):
        '''Similar to ask geno but this function makes sure the translation parameter (the last one) does not go beyond its boundary limits'''
        if number is None or number < 1:
            number = self.sp.popsize
        if xmean is None:
            xmean = self.mean

        if self.countiter == 0:
            self.tic = time.clock()  # backward compatible
            self.elapsed_time = ElapsedTime()

        if self.opts['CMA_AII']:
            if self.countiter == 0:
                self.aii = AII(self.x0, self.sigma0)
            self.flgtelldone = False
            pop = self.aii.ask(number)
            return pop

        sigma = sigma_fac * self.sigma

        # update parameters for sampling the distribution
        #        fac  0      1      10
        # 150-D cigar:
        #           50749  50464   50787
        # 200-D elli:               == 6.9
        #                  99900   101160
        #                 100995   103275 == 2% loss
        # 100-D elli:               == 6.9
        #                 363052   369325  < 2% loss
        #                 365075   365755

        # update distribution
        if self.sp.CMA_on and (
                (self.opts['updatecovwait'] is None and
                 self.countiter >=
                     self.itereigenupdated + 1./(self.sp.c1+self.sp.cmu)/self.N/10
                 ) or
                (self.opts['updatecovwait'] is not None and
                 self.countiter > self.itereigenupdated + self.opts['updatecovwait']
                 )):
            self.updateBD()

        # sample distribution
        if self.flgtelldone:  # could be done in tell()!?
            self.flgtelldone = False
            self.ary = []

        # extract the phenotype and check if the non peridodic decision value is within bounds

        pop_geno = []

        while len(pop_geno) < number:

            # each row is a solution
            arz = self.randn((number, self.N))
    #
            if number == self.sp.popsize:
                self.arz = arz  # is never used
            else:
                pass

            self.ary = self.sigma_vec * np.dot(self.B, (self.D * arz).T).T
            pop = xmean + sigma * self.ary


            for v in pop:
                x = ((v*(space.high-space.low))+space.low )
                
                test = np.any(v < 0)

                for i in xrange(0,len(x),1):
                    if space.high[i] != 360 and not (x[i] < space.low[i] or x[i] > space.high[i]) and len(pop_geno) < number:
                        pop_geno.append( v )

#        print ("POP PHENO WITH REMOVAL "+str(pop_phenoTemp1))

        self.evaluations_per_f_value = 1
        
        # in essence you only return the genotype which phenotype does not exceeds the boundary values for the translation part
        return pop_geno

    def ask(self, number=None, xmean=None, sigma_fac=1):
        ''''''
        pop_geno = self.ask_geno(number, xmean, sigma_fac)



        # N,lambda=20,200: overall CPU 7s vs 5s == 40% overhead, even without bounds!
        #                  new data: 11.5s vs 9.5s == 20%
        # TODO: check here, whether this is necessary?
        # return [self.gp.pheno(x, copy=False, bounds=self.gp.bounds) for x in pop]  # probably fine
        # return [Solution(self.gp.pheno(x, copy=False), copy=False) for x in pop]  # here comes the memory leak, now solved
        # pop_pheno = [Solution(self.gp.pheno(x, copy=False), copy=False).repair(self.gp.bounds) for x in pop_geno]
        pop_pheno = [self.gp.pheno(x, copy=True, bounds=self.gp.bounds) for x in pop_geno]

        if len(pop_pheno) != len(pop_geno):
            print ("pheno-geno not equal")

        if not self.gp.isidentity or use_sent_solutions:  # costs 25% in CPU performance with N,lambda=20,200
            # archive returned solutions, first clean up archive
            if self.countiter % 30/self.popsize**0.5 < 1:
                self.sent_solutions.truncate(0, self.countiter - 1 - 3 * self.N/self.popsize**0.5)
            # insert solutions
            for i in xrange(len(pop_geno)):
                self.sent_solutions[pop_pheno[i]] = {'geno': pop_geno[i],
                                            'pheno': pop_pheno[i],
                                            'iteration': self.countiter}


        return pop_pheno


    # ---------------------------- FUNCTION TO ASK GENO ---------------------------
    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        """"""

        if number is None or number < 1:
            number = self.sp.popsize
        if xmean is None:
            xmean = self.mean

        if self.countiter == 0:
            self.tic = time.clock()  # backward compatible
            self.elapsed_time = ElapsedTime()

        if self.opts['CMA_AII']:
            if self.countiter == 0:
                self.aii = AII(self.x0, self.sigma0)
            self.flgtelldone = False
            pop = self.aii.ask(number)
            return pop

        sigma = sigma_fac * self.sigma

        # update parameters for sampling the distribution
        #        fac  0      1      10
        # 150-D cigar:
        #           50749  50464   50787
        # 200-D elli:               == 6.9
        #                  99900   101160
        #                 100995   103275 == 2% loss
        # 100-D elli:               == 6.9
        #                 363052   369325  < 2% loss
        #                 365075   365755

        # update distribution
        if self.sp.CMA_on and (
                (self.opts['updatecovwait'] is None and
                 self.countiter >=
                     self.itereigenupdated + 1./(self.sp.c1+self.sp.cmu)/self.N/10
                 ) or
                (self.opts['updatecovwait'] is not None and
                 self.countiter > self.itereigenupdated + self.opts['updatecovwait']
                 )):
            self.updateBD()

        # sample distribution
        if self.flgtelldone:  # could be done in tell()!?
            self.flgtelldone = False
            self.ary = []

        # each row is a solution
        arz = self.randn((number, self.N))
        if 11 < 3:  # mutate along the principal axes only
            perm = np.random.permutation(self.N) # indices for mutated principal component
            for i in xrange(min((len(arz), self.N))):
                # perm = np.random.permutation(self.N)  # random principal component, should be much worse
                l = sum(arz[i]**2)**0.5
                arz[i] *= 0
                if 11 < 3: # mirrored sampling
                    arz[i][perm[int(i/2)]] = l * (2 * (i % 2) - 1)
                else:
                    arz[i][perm[i % self.N]] = l * np.sign(np.random.rand(1) - 0.5)
        if number == self.sp.popsize:
            self.arz = arz  # is never used
        else:
            pass

#        if 11 < 3:  # normalize the length to chiN
#            for i in xrange(len(arz)):
#                # arz[i] *= exp(self.randn(1)[0] / 8)
#                ss = sum(arz[i]**2)**0.5
#                arz[i] *= self.const.chiN / ss
#             or to average
#             arz *= 1 * self.const.chiN / np.mean([sum(z**2)**0.5 for z in arz])

        # fac = np.mean(sum(arz**2, 1)**0.5)
        # print fac
        # arz *= self.const.chiN / fac
        self.ary = self.sigma_vec * np.dot(self.B, (self.D * arz).T).T
        pop = xmean + sigma * self.ary
        self.evaluations_per_f_value = 1

        return pop



    # ---------------------------- FUNCTION TO GET MIRROR ---------------------------
    def get_mirror(self, x):
        """"""
        try:
            # dx = x.geno - self.mean, repair or boundary handling is not taken into account
            dx = self.sent_solutions[x]['geno'] - self.mean
        except:
            print('WARNING: use of geno is depreciated')
            dx = self.gp.geno(x, copy=True) - self.mean
        dx *= sum(self.randn(self.N)**2)**0.5 / self.mahalanobisNorm(dx)
        x = self.mean - dx
        y = self.gp.pheno(x, bounds=self.gp.bounds)
        if not self.gp.isidentity or use_sent_solutions:  # costs 25% in CPU performance with N,lambda=20,200
            self.sent_solutions[y] = {'geno': x,
                                        'pheno': y,
                                        'iteration': self.countiter}
        return y

    def mirror_penalized(self, f_values, idx):
        """obsolete and subject to removal (TODO),
        return modified f-values such that for each mirror one becomes worst.

        This function is useless when selective mirroring is applied with no
        more than (lambda-mu)/2 solutions.

        Mirrors are leading and trailing values in ``f_values``.

        """
        assert len(f_values) >= 2 * len(idx)
        m = np.max(np.abs(f_values))
        for i in len(idx):
            if f_values[idx[i]] > f_values[-1-i]:
                f_values[idx[i]] += m
            else:
                f_values[-1-i] += m
        return f_values

    def mirror_idx_cov(self, f_values, idx1):  # will most likely be removed
        ''''''
        idx2 = np.arange(len(f_values) - 1, len(f_values) - 1 - len(idx1), -1)
        f = []
        for i in xrange(len(idx1)):
            f.append(min((f_values[idx1[i]], f_values[idx2[i]])))
            # idx.append(idx1[i] if f_values[idx1[i]] > f_values[idx2[i]] else idx2[i])
        return idx2[np.argsort(f)][-1::-1]

    # ---------------------------- FUNCTION TO ASK AND EVAL ---------------------------
    #
    def ask_and_eval(self, func, args=(), number=None, xmean=None, sigma_fac=1,
                     evaluations=1, aggregation=np.median):
        ''''''
        # initialize
        popsize = self.sp.popsize
        if number is not None:
            popsize = number
        selective_mirroring = True
        nmirrors = self.sp.lam_mirr
        if popsize != self.sp.popsize:
            nmirrors = Mh.sround(popsize * self.sp.lam_mirr / self.sp.popsize)
            # TODO: now selective mirroring might be impaired
        assert nmirrors <= popsize // 2
        self.mirrors_idx = np.arange(nmirrors)  # might never be used
        self.mirrors_rejected_idx = []  # might never be used
        if xmean is None:
            xmean = self.mean

        # do the work
        fit = []  # or np.NaN * np.empty(number)
        X_first = self.ask(popsize)
        X = []
        for k in xrange(int(popsize)):
            nreject = -1
            f = np.NaN
            while f in (np.NaN, None):  # rejection sampling
                nreject += 1
                if k < popsize - nmirrors or nreject:
                    if nreject:
                        x = self.ask(1, xmean, sigma_fac)[0]
                    else:
                        x = X_first.pop(0)
                else:  # mirrored sample
                    if k == popsize - nmirrors and selective_mirroring:
                        self.mirrors_idx = np.argsort(fit)[-1:-1-nmirrors:-1]
                    x = self.get_mirror(X[self.mirrors_idx[popsize - 1 - k]])
                if nreject == 1 and k >= popsize - nmirrors:
                    self.mirrors_rejected_idx.append(k)

                # contraints handling test hardwired ccccccccccc
                if 11 < 3 and self.opts['vv'] and nreject < 2:  # trying out negative C-update as constraints handling
                    if not hasattr(self, 'constraints_paths'):
                        k = 1
                        self.constraints_paths = [np.zeros(self.N) for _i in xrange(k)]
                    Izero = np.zeros([self.N, self.N])
                    for i in xrange(self.N):
                        if x[i] < 0:
                            Izero[i][i] = 1
                            self.C -= self.opts['vv'] * Izero
                            Izero[i][i] = 0
                    if 1 < 3 and sum([ (9 + i + 1) * x[i] for i in xrange(self.N)]) > 50e3:
                        self.constraints_paths[0] = 0.9 * self.constraints_paths[0] + 0.1 * (x - self.mean) / self.sigma
                        self.C -= (self.opts['vv'] / self.N) * np.outer(self.constraints_paths[0], self.constraints_paths[0])

                f = func(x, *args)
                if f not in (np.NaN, None) and evaluations > 1:
                    f = aggregation([f] + [func(x, *args) for _i in xrange(int(evaluations-1))])
                if nreject + 1 % 1000 == 0:
                    print('  %d solutions rejected (f-value NaN or None) at iteration %d' %
                          (nreject, self.countiter))
            fit.append(f)
            X.append(x)
        self.evaluations_per_f_value = int(evaluations)
        return X, fit

    def tell(self, solutions, function_values, check_points=None, copy=False):
        """pass objective function values to prepare for next
        iteration. This core procedure of the CMA-ES algorithm updates
        all state variables, in particular the two evolution paths, the
        distribution mean, the covariance matrix and a step-size.

        Arguments
        ---------
            `solutions`
                list or array of candidate solution points (of
                type `numpy.ndarray`), most presumably before
                delivered by method `ask()` or `ask_and_eval()`.
            `function_values`
                list or array of objective function values
                corresponding to the respective points. Beside for termination
                decisions, only the ranking of values in `function_values`
                is used.
            `check_points`
                If ``check_points is None``, only solutions that are not generated
                by `ask()` are possibly clipped (recommended). ``False`` does not clip
                any solution (not recommended).
                If ``True``, clips solutions that realize long steps (i.e. also
                those that are unlikely to be generated with `ask()`). `check_points`
                can be a list of indices to be checked in solutions.
            `copy`
                ``solutions`` can be modified in this routine, if ``copy is False``

        Details
        -------
        `tell()` updates the parameters of the multivariate
        normal search distribution, namely covariance matrix and
        step-size and updates also the attributes `countiter` and
        `countevals`. To check the points for consistency is quadratic
        in the dimension (like sampling points).

        Bugs
        ----
        The effect of changing the solutions delivered by `ask()` depends on whether
        boundary handling is applied. With boundary handling, modifications are
        disregarded. This is necessary to apply the default boundary handling that
        uses unrepaired solutions but might change in future.

        Example
        -------
        ::

            import cma
            func = cma.fcts.elli  # choose objective function
            es = cma.CMAEvolutionStrategy(cma.np.random.rand(10), 1)
            while not es.stop():
               X = es.ask()
               es.tell(X, [func(x) for x in X])
            es.result()  # where the result can be found

        :See: class `CMAEvolutionStrategy`, `ask()`, `ask_and_eval()`, `fmin()`

        """
    #____________________________________________________________
    # TODO: consider an input argument that flags injected trust-worthy solutions (which means
    #       that they can be treated "absolut" rather than "relative")

        self.unsuccess = 0

        if self.flgtelldone:
            raise _Error('tell should only be called once per iteration')

        lam = len(solutions)
        if lam != array(function_values).shape[0]:
            raise _Error('for each candidate solution '
                        + 'a function value must be provided')
        if lam + self.sp.lam_mirr < 3:
            raise _Error('population size ' + str(lam) + ' is too small when option CMA_mirrors * popsize < 0.5')

        if not np.isscalar(function_values[0]):
            if np.isscalar(function_values[0][0]):
                if self.countiter <= 1:
                    print('WARNING: function values are not a list of scalars (further warnings are suppressed)')
                function_values = [val[0] for val in function_values]
            else:
                raise _Error('objective function values must be a list of scalars')


        ### prepare
        N = self.N
        sp = self.sp
        if 11 < 3 and lam != sp.popsize:  # turned off, because mu should stay constant, still not desastrous
            print('WARNING: population size has changed, recomputing parameters')
            self.sp.set(self.opts, lam)  # not really tested
        if lam < sp.mu:  # rather decrease cmean instead of having mu > lambda//2
            raise _Error('not enough solutions passed to function tell (mu>lambda)')

        self.countiter += 1  # >= 1 now
        self.countevals += sp.popsize * self.evaluations_per_f_value
        self.best.update(solutions, self.sent_solutions, function_values, self.countevals)

        flgseparable = self.opts['CMA_diagonal'] is True \
                       or self.countiter <= self.opts['CMA_diagonal']
        if not flgseparable and len(self.C.shape) == 1:  # C was diagonal ie 1-D
            # enter non-separable phase (no easy return from here)
            self.B = np.eye(N) # identity(N)
            self.C = np.diag(self.C)
            idx = np.argsort(self.D)
            self.D = self.D[idx]
            self.B = self.B[:,idx]
            self.Zneg = np.zeros((N, N))

        ### manage fitness
        fit = self.fit  # make short cut

        # CPU for N,lam=20,200: this takes 10s vs 7s
        fit.bndpen = self.boundPenalty.update(function_values, self)(solutions, self.sent_solutions, self.gp)
        # for testing:
        # fit.bndpen = self.boundPenalty.update(function_values, self)([s.unrepaired for s in solutions])
        fit.idx = np.argsort(array(fit.bndpen) + array(function_values))
        fit.fit = array(function_values, copy=False)[fit.idx]

        # update output data TODO: this is obsolete!? However: need communicate current best x-value?
        # old: out['recent_x'] = self.gp.pheno(pop[0])
        self.out['recent_x'] = array(solutions[fit.idx[0]])  # TODO: change in a data structure(?) and use current as identify
        self.out['recent_f'] = fit.fit[0]

        # fitness histories
        fit.hist.insert(0, fit.fit[0])
        # if len(self.fit.histbest) < 120+30*N/sp.popsize or  # does not help, as tablet in the beginning is the critical counter-case
        if ((self.countiter % 5) == 0):  # 20 percent of 1e5 gen.
            fit.histbest.insert(0, fit.fit[0])
            fit.histmedian.insert(0, np.median(fit.fit) if len(fit.fit) < 21
                                    else fit.fit[self.popsize // 2])
        if len(fit.histbest) > 2e4: # 10 + 30*N/sp.popsize:
            fit.histbest.pop()
            fit.histmedian.pop()
        if len(fit.hist) > 10 + 30*N/sp.popsize:
            fit.hist.pop()

        if self.opts['CMA_AII']:
            self.aii.tell(solutions, function_values)
            self.flgtelldone = True
            # for output:
            self.mean = self.aii.mean
            self.dC = self.aii.sigmai**2
            self.sigma = self.aii.sigma
            self.D = 1e-11 + (self.aii.r**2)**0.5
            self.more_to_write += [self.aii.sigma_r]
            return

        # TODO: clean up inconsistency when an unrepaired solution is available and used
        pop = []  # create pop from input argument solutions
        for s in solutions:  # use phenotype before Solution.repair()
            if use_sent_solutions:
                x = self.sent_solutions.pop(s, None)  # 12.7s vs 11.3s with N,lambda=20,200
                if x is not None:
                    pop.append(x['geno'])
                    # TODO: keep additional infos or don't pop s from sent_solutions in the first place
                else:
                    # print 'WARNING: solution not found in ``self.sent_solutions`` (is expected for injected solutions)'
                    pop.append(self.gp.geno(s, copy=copy))  # cannot recover the original genotype with boundary handling
                    if check_points in (None, True, 1):
                        self.repair_genotype(pop[-1])  # necessary if pop[-1] was changed or injected by the user.
            else:  # TODO: to be removed?
                # print 'WARNING: ``geno`` mapping depreciated'
                pop.append(self.gp.geno(s, copy=copy))
                if check_points in (None, True, 1):
                    self.repair_genotype(pop[-1])  # necessary or not?
                # print 'repaired'

        mold = self.mean
        sigma_fac = 1

        # check and normalize each x - m
        # check_points is a flag (None is default: check non-known solutions) or an index list
        # should also a number possible (first check_points points)?
        if check_points not in (None, False, 0, [], ()):  # useful in case of injected solutions and/or adaptive encoding, however is automatic with use_sent_solutions
            try:
                if len(check_points):
                    idx = check_points
            except:
                idx = xrange(sp.popsize)

            for k in idx:
                self.repair_genotype(pop[k])

        # sort pop
        if type(pop) is not array: # only arrays can be multiple indexed
            pop = array(pop, copy=False)

        pop = pop[fit.idx]

        if self.opts['CMA_elitist'] and self.best.f < fit.fit[0]:
            if self.best.x_geno is not None:
                xp = [self.best.x_geno]
                # xp = [self.best.xdict['geno']]
                # xp = [self.gp.geno(self.best.x[:])]  # TODO: remove
                # print self.mahalanobisNorm(xp[0]-self.mean)
                self.clip_or_fit_solutions(xp, [0])
                pop = array([xp[0]] + list(pop))
            else:
                print('genotype for elitist not found')

        # compute new mean
        self.mean = mold + self.sp.cmean * \
                    (sum(sp.weights * pop[0:sp.mu].T, 1) - mold)


        # check Delta m (this is not default, but could become at some point)
        # CAVE: upper_length=sqrt(2)+2 is too restrictive, test upper_length = sqrt(2*N) thoroughly.
        # simple test case injecting self.mean:
        # self.mean = 1e-4 * self.sigma * np.random.randn(N)
        if 11 < 3 and self.opts['vv'] and check_points:  # TODO: check_points might be an index-list
            cmean = self.sp.cmean / min(1, (sqrt(self.opts['vv']*N)+2) / ( # abuse of cmean
                (sqrt(self.sp.mueff) / self.sp.cmean) *
                self.mahalanobisNorm(self.mean - mold)))
        else:
            cmean = self.sp.cmean

        if 11 < 3:  # plot length of mean - mold
            self.more_to_write += [sqrt(sp.mueff) *
                sum(((1./self.D) * dot(self.B.T, self.mean - mold))**2)**0.5 /
                       self.sigma / sqrt(N) / cmean]

        # get learning rate constants
        cc, c1, cmu = sp.cc, sp.c1, sp.cmu
        if flgseparable:
            cc, c1, cmu = sp.cc_sep, sp.c1_sep, sp.cmu_sep

        # now the real work can start

        # evolution paths
        self.ps = (1-sp.cs) * self.ps + \
                  (sqrt(sp.cs*(2-sp.cs)*sp.mueff)  / self.sigma / cmean) * \
                  dot(self.B, (1./self.D) * dot(self.B.T, (self.mean - mold) / self.sigma_vec))

        # "hsig", correction with self.countiter seems not necessary, also pc starts with zero
        hsig = sum(self.ps**2) / (1-(1-sp.cs)**(2*self.countiter)) / self.N < 2 + 4./(N+1)
        if 11 < 3:
            # hsig = 1
            # sp.cc = 4 / (N + 4)
            # sp.cs = 4 / (N + 4)
            # sp.cc = 1
            # sp.damps = 2  #
            # sp.CMA_on = False
            # c1 = 0  # 2 / ((N + 1.3)**2 + 0 * sp.mu) # 1 / N**2
            # cmu = min([1 - c1, cmu])
            if self.countiter == 1:
                print('parameters modified')
        # hsig = sum(self.ps**2) / self.N < 2 + 4./(N+1)
        # adjust missing variance due to hsig, in 4-D with damps=1e99 and sig0 small
        #       hsig leads to premature convergence of C otherwise
        #hsiga = (1-hsig**2) * c1 * cc * (2-cc)  # to be removed in future
        c1a = c1 - (1-hsig**2) * c1 * cc * (2-cc)  # adjust for variance loss

        if 11 < 3:  # diagnostic data
            self.out['hsigcount'] += 1 - hsig
            if not hsig:
                self.hsiglist.append(self.countiter)
        if 11 < 3:  # diagnostic message
            if not hsig:
                print(str(self.countiter) + ': hsig-stall')
        if 11 < 3:  # for testing purpose
            hsig = 1 # TODO:
            #       put correction term, but how?
            if self.countiter == 1:
                print('hsig=1')

        self.pc = (1-cc) * self.pc + \
                  hsig * (sqrt(cc*(2-cc)*sp.mueff) / self.sigma / cmean) * \
                  (self.mean - mold)  / self.sigma_vec

        # covariance matrix adaptation/udpate
        if sp.CMA_on:
            # assert sp.c1 + sp.cmu < sp.mueff / N  # ??
            assert c1 + cmu <= 1

            # default full matrix case
            if not flgseparable:
                Z = (pop[0:sp.mu] - mold) / (self.sigma * self.sigma_vec)
                Z = dot((cmu * sp.weights) * Z.T, Z)  # learning rate integrated
                if self.sp.neg.cmuexp:
                    tmp = (pop[-sp.neg.mu:] - mold) / (self.sigma * self.sigma_vec)
                    self.Zneg *= 1 - self.sp.neg.cmuexp  # for some reason necessary?
                    self.Zneg += dot(sp.neg.weights * tmp.T, tmp) - self.C
                    # self.update_exponential(dot(sp.neg.weights * tmp.T, tmp) - 1 * self.C, -1*self.sp.neg.cmuexp)

                if 11 < 3: # ?3 to 5 times slower??
                    Z = np.zeros((N,N))
                    for k in xrange(sp.mu):
                        z = (pop[k]-mold)
                        Z += np.outer((cmu * sp.weights[k] / (self.sigma * self.sigma_vec)**2) * z, z)

                self.C *= 1 - c1a - cmu
                self.C += np.outer(c1 * self.pc, self.pc) + Z
                self.dC = np.diag(self.C)  # for output and termination checking

            else: # separable/diagonal linear case
                assert(c1+cmu <= 1)
                Z = np.zeros(N)
                for k in xrange(sp.mu):
                    z = (pop[k]-mold) / (self.sigma * self.sigma_vec) # TODO see above
                    Z += sp.weights[k] * z * z  # is 1-D
                self.C = (1-c1a-cmu) * self.C + c1 * self.pc * self.pc + cmu * Z
                # TODO: self.C *= exp(cmuneg * (N - dot(sp.neg.weights,  **2)
                self.dC = self.C
                self.D = sqrt(self.C)  # C is a 1-D array
                self.itereigenupdated = self.countiter

                # idx = self.mirror_idx_cov()  # take half of mirrored vectors for negative update

        # qqqqqqqqqqq
        if 1 < 3 and np.isfinite(sp.dampsvec):
            if self.countiter == 1:
                print("WARNING: CMA_dampsvec option is experimental")
            sp.dampsvec *= np.exp(sp.dampsvec_fading/self.N)
            # TODO: rank-lambda update: *= (1 + sum(z[z>1]**2-1) * exp(sum(z[z<1]**2-1))
            self.sigma_vec *= np.exp((sp.cs/sp.dampsvec/2) * (self.ps**2 - 1))
            # self.sigma_vec *= np.exp((sp.cs/sp.dampsvec) * (abs(self.ps) - (2/np.pi)**0.5))
            self.more_to_write += [exp(np.mean((self.ps**2 - 1)**2))]
            # TODO: rank-mu update

        # step-size adaptation, adapt sigma
        if 1 < 3:  #
            self.sigma *= sigma_fac * \
                            np.exp((min((1, (sp.cs/sp.damps) *
                                    (sqrt(sum(self.ps**2))/self.const.chiN - 1)))))
        else:
            self.sigma *= sigma_fac * \
                            np.exp((min((1000, (sp.cs/sp.damps/2) *
                                    (sum(self.ps**2)/N - 1)))))
        if 11 < 3:
            # derandomized MSR = natural gradient descent using mean(z**2) instead of mu*mean(z)**2
            lengths = array([sum(z**2)**0.5 for z in self.arz[fit.idx[:self.sp.mu]]])
            # print lengths[0::int(self.sp.mu/5)]
            self.sigma *= np.exp(self.sp.mueff**0.5 * dot(self.sp.weights, lengths / self.const.chiN - 1))**(2/(N+1))

        if 11 < 3 and self.opts['vv']:
            if self.countiter < 2:
                print('constant sigma applied')
                print(self.opts['vv'])  # N=10,lam=10: 0.8 is optimal
            self.sigma = self.opts['vv'] * self.sp.mueff * sum(self.mean**2)**0.5 / N

        if self.sigma * min(self.dC)**0.5 < self.opts['minstd']:
            self.sigma = self.opts['minstd'] / min(self.dC)**0.5
        # g = self.countiter
        # N = self.N
        mindx = eval(self.opts['mindx']) if type(self.opts['mindx']) == type('') else self.opts['mindx']
        if self.sigma * min(self.D) < mindx:  # TODO: sigma_vec is missing here
            self.sigma = mindx / min(self.D)

        if self.sigma > 1e9 * self.sigma0:
            alpha = self.sigma / max(self.D)
            self.multiplyC(alpha)
            self.sigma /= alpha**0.5
            self.opts['tolupsigma'] /= alpha**0.5  # to be compared with sigma

        # TODO increase sigma in case of a plateau?

        # Uncertainty noise measurement is done on an upper level

        # output, has moved up, e.g. as part of fmin, TODO to be removed
        if 11 < 3 and self.opts['verb_log'] > 0 and (self.countiter < 4 or
                                          self.countiter % self.opts['verb_log'] == 0):
            # this assumes that two logger with the same name access the same data!
            CMADataLogger(self.opts['verb_filenameprefix']).register(self, append=True).add()
            # self.writeOutput(solutions[fit.idx[0]])

        self.flgtelldone = True
    # end tell()

    #____________________________________________________________
    def tell_ViE(self, solutions, function_values, suvivorsNo, unsuccess, check_points=None, copy=False):
        ''''''
    #____________________________________________________________
    # TODO: consider an input argument that flags injected trust-worthy solutions (which means
    #       that they can be treated "absolut" rather than "relative")

        self.unsuccess = unsuccess * 12


        if self.flgtelldone:
            raise _Error('tell should only be called once per iteration')

        lam = len(solutions)
        if lam != array(function_values).shape[0]:
            raise _Error('for each candidate solution '
                        + 'a function value must be provided')
#        if lam + self.sp.lam_mirr < 3:
#            raise _Error('population size ' + str(lam) + ' is too small when option CMA_mirrors * popsize < 0.5')

        if not np.isscalar(function_values[0]):
            if np.isscalar(function_values[0][0]):
                if self.countiter <= 1:
                    print('WARNING: function values are not a list of scalars (further warnings are suppressed)')
                function_values = [val[0] for val in function_values]
            else:
                raise _Error('objective function values must be a list of scalars')

        ### prepare
        N = self.N # the number of dimensions
        if 1 < 3:
            self.sp = CMAViEParameters(N, self.opts, suvivorsNo)
        sp = self.sp # the theoretical popsize


#        if 11 < 3 and lam != sp.popsize:  # turned off, because mu should stay constant, still not desastrous
#            print('WARNING: population size has changed, recomputing parameters')
#            self.sp.set(self.opts, lam)  # not really tested

        # ------------------ CASE WHEN LAMBDA IS LOWER THAN THE MINIMUM REQUIRED TO COMPUTE MEAN
#        if lam < sp.mu:  # rather decrease cmean instead of having mu > lambda//2
#            print ("lambda: "+str(lam))
#            print ("mu: "+str(sp.mu))
#            raise _Error('not enough solutions passed to function tell (mu>lambda)')

        # ---------------------- COUNTING ITERATIONS AND EVALSself.unsucess
        self.countiter += 1  # >= 1 now
        self.countevals += sp.popsize * self.evaluations_per_f_value
        self.best.update(solutions, self.sent_solutions, function_values, self.countevals)

        flgseparable = self.opts['CMA_diagonal'] is True \
                       or self.countiter <= self.opts['CMA_diagonal']
        if not flgseparable and len(self.C.shape) == 1:  # C was diagonal ie 1-D
            # enter non-separable phase (no easy return from here)
            self.B = np.eye(N) # identity(N)
            self.C = np.diag(self.C)
            idx = np.argsort(self.D)
            self.D = self.D[idx]
            self.B = self.B[:,idx]
            self.Zneg = np.zeros((N, N))

        ### manage fitness
        fit = self.fit  # make short cut, this fit contains the fitness values of the previous generations

        # CPU for N,lam=20,200: this takes 10s vs 7s

        # ----------------- UPDATING BOUNDARY PENALTIES
        fit.bndpen = self.boundPenalty.update(function_values, self)(solutions, self.sent_solutions, self.gp)


         # ---------------- SORTING THE NEW FITNESS VALUES
        fit.idx = np.argsort(array(fit.bndpen) + array(function_values))
        fit.fit = array(function_values, copy=False)[fit.idx]


        # ----------------- UPDATE THE LATEST BEST SOLUTION
        self.out['recent_x'] = array(solutions[fit.idx[0]])  # TODO: change in a data structure(?) and use current as identify
        self.out['recent_f'] = fit.fit[0]


        # ----------------- UPDATE FITNESS HISTORY
        fit.hist.insert(0, fit.fit[0])
        # if len(self.fit.histbest) < 120+30*N/sp.popsize or  # does not help, as tablet in the beginning is the critical counter-case
        if ((self.countiter % 5) == 0):  # 20 percent of 1e5 gen.
            fit.histbest.insert(0, fit.fit[0])
            fit.histmedian.insert(0, np.median(fit.fit) if len(fit.fit) < 21
                                    else fit.fit[self.popsize // 2])
        if len(fit.histbest) > 2e4: # 10 + 30*N/sp.popsize:
            fit.histbest.pop()
            fit.histmedian.pop()
        if len(fit.hist) > 10 + 30*N/sp.popsize:
            fit.hist.pop()


        # --------------------- this one is not used
        if self.opts['CMA_AII']:
            self.aii.tell(solutions, function_values)
            self.flgtelldone = True
            # for output:
            self.mean = self.aii.mean
            self.dC = self.aii.sigmai**2
            self.sigma = self.aii.sigma
            self.D = 1e-11 + (self.aii.r**2)**0.5
            self.more_to_write += [self.aii.sigma_r]
            return

        # --------------------- LOADING THE SOLUTIONS GENOTYPES INTO pop array
        # TODO: clean up inconsistency when an unrepaired solution is available and used
        pop = []  # create pop from input argument solutions || contains the genotypes of the solutions
        for s in solutions:  # use phenotype before Solution.repair()
            if use_sent_solutions: # this is a boolean that seems to be always True
                x = self.sent_solutions.pop(s, None)  # 12.7s vs 11.3s with N,lambda=20,200 || x contains the hash geno, pheno, iterationno
                if x is not None:
                    pop.append(x['geno'])
                    # TODO: keep additional infos or don't pop s from sent_solutions in the first place
                else:
                    # print 'WARNING: solution not found in ``self.sent_solutions`` (is expected for injected solutions)'
                    pop.append(self.gp.geno(s, copy=copy))  # cannot recover the original genotype with boundary handling
                    if check_points in (None, True, 1):
                        self.repair_genotype(pop[-1])  # necessary if pop[-1] was changed or injected by the user.
#            else:  # TODO: to be removed?
#                # print 'WARNING: ``geno`` mapping depreciated'
#                pop.append(self.gp.geno(s, copy=copy))
#                if check_points in (None, True, 1):
#                    self.repair_genotype(pop[-1])  # necessary or not?
#                # print 'repaired'

        mold = self.mean # simply the mean -> [0.1 0.1 0.1 0.1 0.1] for 5D for example
        sigma_fac = 1

        # --------------------- not going in this so far
        # check and normalize each x - m
        # check_points is a flag (None is default: check non-known solutions) or an index list
        # should also a number possible (first check_points points)?
        if check_points not in (None, False, 0, [], ()):  # useful in case of injected solutions and/or adaptive encoding, however is automatic with use_sent_solutions
            try:
                if len(check_points):
                    idx = check_points
            except:
                idx = xrange(sp.popsize)

            for k in idx:
                self.repair_genotype(pop[k])

        # --------------------- SORTING THE REST OF THE SOLUTIONS BY THEIR FITNESS VALUES (the fit.idx found before)
        if type(pop) is not array: # only arrays can be multiple indexed
            pop = array(pop, copy=False)

        pop = pop[fit.idx]

        # --------------------- the CMA_elitist case, not of concern for now
        if self.opts['CMA_elitist'] and self.best.f < fit.fit[0]:
            if self.best.x_geno is not None:
                xp = [self.best.x_geno]
                # xp = [self.best.xdict['geno']]
                # xp = [self.gp.geno(self.best.x[:])]  # TODO: remove
                # print self.mahalanobisNorm(xp[0]-self.mean)
                self.clip_or_fit_solutions(xp, [0])
                pop = array([xp[0]] + list(pop))
            else:
                print('genotype for elitist not found')

        # --------------------- COMPUTING THE NEW MEAN
        # because the old mean was calculated based on the sp.mu which is fixed we need to remove it
        # self.sp.mean is alway 1.0, is it a scaling factor of some sort?

        #self.mean = mold + self.sp.cmean * (sum(sp.weights * pop[0:sp.mu].T, 1) - mold)

        if sp.wholePopFlag:
            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop.T, 1) - mold)
        else:
            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop[0:sp.mu].T, 1) - mold)
#            self.mean = self.sp.cmean * sum(sp.weights * pop[0:sp.mu].T, 1) # this mean is modified with the +- mold
#            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop.T, 1) - mold)
#        print(">> new mean: "+str(self.mean))


        # check Delta m (this is not default, but could become at some point)
        # CAVE: upper_length=sqrt(2)+2 is too restrictive, test upper_length = sqrt(2*N) thoroughly.
        # simple test case injecting self.mean:
        # self.mean = 1e-4 * self.sigma * np.random.randn(N)
        if 11 < 3 and self.opts['vv'] and check_points:  # TODO: check_points might be an index-list
            cmean = self.sp.cmean / min(1, (sqrt(self.opts['vv']*N)+2) / ( # abuse of cmean
                (sqrt(self.sp.mueff) / self.sp.cmean) *
                self.mahalanobisNorm(self.mean - mold)))
        else:
            cmean = self.sp.cmean

        if 11 < 3:  # plot length of mean - mold
            self.more_to_write += [sqrt(sp.mueff) *
                sum(((1./self.D) * dot(self.B.T, self.mean - mold))**2)**0.5 /
                       self.sigma / sqrt(N) / cmean]

        # get learning rate constants
        cc, c1, cmu = sp.cc, sp.c1, sp.cmu
        if flgseparable:
            cc, c1, cmu = sp.cc_sep, sp.c1_sep, sp.cmu_sep

        # now the real work can start

        # evolution paths
        self.ps = (1-sp.cs) * self.ps + \
                  (sqrt(sp.cs*(2-sp.cs)*sp.mueff)  / self.sigma / cmean) * \
                  dot(self.B, (1./self.D) * dot(self.B.T, (self.mean - mold) / self.sigma_vec))

        # "hsig", correction with self.countiter seems not necessary, also pc starts with zero
        hsig = sum(self.ps**2) / (1-(1-sp.cs)**(2*self.countiter)) / self.N < 2 + 4./(N+1)
        if 11 < 3:
            # hsig = 1
            # sp.cc = 4 / (N + 4)
            # sp.cs = 4 / (N + 4)
            # sp.cc = 1
            # sp.damps = 2  #
            # sp.CMA_on = False
            # c1 = 0  # 2 / ((N + 1.3)**2 + 0 * sp.mu) # 1 / N**2
            # cmu = min([1 - c1, cmu])
            if self.countiter == 1:
                print('parameters modified')
        # hsig = sum(self.ps**2) / self.N < 2 + 4./(N+1)
        # adjust missing variance due to hsig, in 4-D with damps=1e99 and sig0 small
        #       hsig leads to premature convergence of C otherwise
        #hsiga = (1-hsig**2) * c1 * cc * (2-cc)  # to be removed in future
        c1a = c1 - (1-hsig**2) * c1 * cc * (2-cc)  # adjust for variance loss

        if 11 < 3:  # diagnostic data
            self.out['hsigcount'] += 1 - hsig
            if not hsig:
                self.hsiglist.append(self.countiter)
        if 11 < 3:  # diagnostic message
            if not hsig:
                print(str(self.countiter) + ': hsig-stall')
        if 11 < 3:  # for testing purpose
            hsig = 1 # TODO:
            #       put correction term, but how?
            if self.countiter == 1:
                print('hsig=1')

        self.pc = (1-cc) * self.pc + \
                  hsig * (sqrt(cc*(2-cc)*sp.mueff) / self.sigma / cmean) * \
                  (self.mean - mold)  / self.sigma_vec

        # covariance matrix adaptation/udpate
        if sp.CMA_on:
            # assert sp.c1 + sp.cmu < sp.mueff / N  # ??
            assert c1 + cmu <= 1

            # default full matrix case
            if not flgseparable:
                if sp.wholePopFlag:
                    Z = (pop - mold) / (self.sigma * self.sigma_vec)
                else:
                    Z = (pop[0:sp.mu] - mold) / (self.sigma * self.sigma_vec)
                Z = dot((cmu * sp.weights) * Z.T, Z)  # learning rate integrated
                if self.sp.neg.cmuexp:
                    tmp = (pop[-sp.neg.mu:] - mold) / (self.sigma * self.sigma_vec)
                    self.Zneg *= 1 - self.sp.neg.cmuexp  # for some reason necessary?
                    self.Zneg += dot(sp.neg.weights * tmp.T, tmp) - self.C
                    # self.update_exponential(dot(sp.neg.weights * tmp.T, tmp) - 1 * self.C, -1*self.sp.neg.cmuexp)

#                if 11 < 3: # ?3 to 5 times slower??
#                    Z = np.zeros((N,N))
#                    for k in xrange(sp.mu):
#                        z = (pop[k]-mold)
#                        Z += np.outer((cmu * sp.weights[k] / (self.sigma * self.sigma_vec)**2) * z, z)

                self.C *= 1 - c1a - cmu
                self.C += np.outer(c1 * self.pc, self.pc) + Z
#                self.C = np.eye(self.N, self.N);    # NEW , remove it later !!!!!!!! ILYA
                self.dC = np.diag(self.C)  # for output and termination checking

            else: # separable/diagonal linear case
                assert(c1+cmu <= 1)
                Z = np.zeros(N)
                for k in xrange(sp.mu):
                    z = (pop[k]-mold) / (self.sigma * self.sigma_vec) # TODO see above
                    Z += sp.weights[k] * z * z  # is 1-D
                self.C = (1-c1a-cmu) * self.C + c1 * self.pc * self.pc + cmu * Z
                # TODO: self.C *= exp(cmuneg * (N - dot(sp.neg.weights,  **2)
#                self.C = np.eye(self.N, self.N);   # NEW , remove it later !!!!!!!! ILYA
                self.dC = self.C
                self.D = sqrt(self.C)  # C is a 1-D array
                self.itereigenupdated = self.countiter

                # idx = self.mirror_idx_cov()  # take half of mirrored vectors for negative update

        # qqqqqqqqqqq
        if 1 < 3 and np.isfinite(sp.dampsvec):
            if self.countiter == 1:
                print("WARNING: CMA_dampsvec option is experimental")
            sp.dampsvec *= np.exp(sp.dampsvec_fading/self.N)
            # TODO: rank-lambda update: *= (1 + sum(z[z>1]**2-1) * exp(sum(z[z<1]**2-1))
            self.sigma_vec *= np.exp((sp.cs/sp.dampsvec/2) * (self.ps**2 - 1))
            # self.sigma_vec *= np.exp((sp.cs/sp.dampsvec) * (abs(self.ps) - (2/np.pi)**0.5))
            self.more_to_write += [exp(np.mean((self.ps**2 - 1)**2))]
            # TODO: rank-mu update

        # step-size adaptation, adapt sigma
        if 1 < 3:  #
            self.sigma *= sigma_fac * \
                            np.exp((min((1, (sp.cs/sp.damps) *
                                    (sqrt(sum(self.ps**2))/self.const.chiN - 1)))))
        else:
            self.sigma *= sigma_fac * \
                            np.exp((min((1000, (sp.cs/sp.damps/2) *
                                    (sum(self.ps**2)/N - 1)))))
#        if 11 < 3:
#            # derandomized MSR = natural gradient descent using mean(z**2) instead of mu*mean(z)**2
#            lengths = array([sum(z**2)**0.5 for z in self.arz[fit.idx[:self.sp.mu]]])
#            # print lengths[0::int(self.sp.mu/5)]
#            self.sigma *= np.exp(self.sp.mueff**0.5 * dot(self.sp.weights, lengths / self.const.chiN - 1))**(2/(N+1))
#
#        if 11 < 3 and self.opts['vv']:
#            if self.countiter < 2:
#                print('constant sigma applied')
#                print(self.opts['vv'])  # N=10,lam=10: 0.8 is optimal
#            self.sigma = self.opts['vv'] * self.sp.mueff * sum(self.mean**2)**0.5 / N

        if self.sigma * min(self.dC)**0.5 < self.opts['minstd']:
            self.sigma = self.opts['minstd'] / min(self.dC)**0.5
        # g = self.countiter
        # N = self.N
        mindx = eval(self.opts['mindx']) if type(self.opts['mindx']) == type('') else self.opts['mindx']
        if self.sigma * min(self.D) < mindx:  # TODO: sigma_vec is missing here
            self.sigma = mindx / min(self.D)

        if self.sigma > 1e9 * self.sigma0:
            alpha = self.sigma / max(self.D)
            self.multiplyC(alpha)
            self.sigma /= alpha**0.5
            self.opts['tolupsigma'] /= alpha**0.5  # to be compared with sigma

        # TODO increase sigma in case of a plateau?

        # Uncertainty noise measurement is done on an upper level

#        # output, has moved up, e.g. as part of fmin, TODO to be removed
#        if 11 < 3 and self.opts['verb_log'] > 0 and (self.countiter < 4 or
#                                          self.countiter % self.opts['verb_log'] == 0):
#            # this assumes that two logger with the same name access the same data!
#            CMADataLogger(self.opts['verb_filenameprefix']).register(self, append=True).add()
#            # self.writeOutput(solutions[fit.idx[0]])

        self.flgtelldone = True
    # end tell()

    def tell_ViE_Constraint(self, solutions, function_values, indeces, unsuccess, check_points=None, copy=False):
        ''''''
    #____________________________________________________________
    # TODO: consider an input argument that flags injected trust-worthy solutions (which means
    #       that they can be treated "absolut" rather than "relative")




        if self.flgtelldone:
            raise _Error('tell should only be called once per iteration')

        lam = len(solutions)
        if lam != array(function_values).shape[0]:
            raise _Error('for each candidate solution '
                        + 'a function value must be provided')
#        if lam + self.sp.lam_mirr < 3:
#            raise _Error('population size ' + str(lam) + ' is too small when option CMA_mirrors * popsize < 0.5')

        if not np.isscalar(function_values[0]):
            if np.isscalar(function_values[0][0]):
                if self.countiter <= 1:
                    print('WARNING: function values are not a list of scalars (further warnings are suppressed)')
                function_values = [val[0] for val in function_values]
            else:
                raise _Error('objective function values must be a list of scalars')

        ### prepare
        N = self.N # the number of dimensions
        if 1 < 3:
            self.sp = CMAViEConsParameters(N, self.opts, len(indeces[1]))
        sp = self.sp # the theoretical popsize

        self.unsuccess = unsuccess * sp.popsize


#        if 11 < 3 and lam != sp.popsize:  # turned off, because mu should stay constant, still not desastrous
#            print('WARNING: population size has changed, recomputing parameters')
#            self.sp.set(self.opts, lam)  # not really tested

        # ------------------ CASE WHEN LAMBDA IS LOWER THAN THE MINIMUM REQUIRED TO COMPUTE MEAN
#        if lam < sp.mu:  # rather decrease cmean instead of having mu > lambda//2
#            print ("lambda: "+str(lam))
#            print ("mu: "+str(sp.mu))
#            raise _Error('not enough solutions passed to function tell (mu>lambda)')

        # ---------------------- COUNTING ITERATIONS AND EVALSself.unsucess
        self.countiter += 1  # >= 1 now
        self.countevals += sp.popsize * self.evaluations_per_f_value
        self.best.update(solutions, self.sent_solutions, function_values, self.countevals)

        flgseparable = self.opts['CMA_diagonal'] is True \
                       or self.countiter <= self.opts['CMA_diagonal']
        if not flgseparable and len(self.C.shape) == 1:  # C was diagonal ie 1-D
            # enter non-separable phase (no easy return from here)
            self.B = np.eye(N) # identity(N)
            self.C = np.diag(self.C)
            idx = np.argsort(self.D)
            self.D = self.D[idx]
            self.B = self.B[:,idx]
            self.Zneg = np.zeros((N, N))

        ### manage fitness
        fit = self.fit  # make short cut, this fit contains the fitness values of the previous generations

        # CPU for N,lam=20,200: this takes 10s vs 7s

        # ----------------- UPDATING BOUNDARY PENALTIES
        fit.bndpen = self.boundPenalty.update(function_values, self)(solutions, self.sent_solutions, self.gp)


         # ---------------- SORTING THE NEW FITNESS VALUES
        fit.idx = np.argsort(array(fit.bndpen) + array(function_values))
        fit.fit = array(function_values, copy=False)[fit.idx]


        # ----------------- UPDATE THE LATEST BEST SOLUTION
        self.out['recent_x'] = array(solutions[fit.idx[0]])  # TODO: change in a data structure(?) and use current as identify
        self.out['recent_f'] = fit.fit[0]


        # ----------------- UPDATE FITNESS HISTORY
        fit.hist.insert(0, fit.fit[0])
        # if len(self.fit.histbest) < 120+30*N/sp.popsize or  # does not help, as tablet in the beginning is the critical counter-case
        if ((self.countiter % 5) == 0):  # 20 percent of 1e5 gen.
            fit.histbest.insert(0, fit.fit[0])
            fit.histmedian.insert(0, np.median(fit.fit) if len(fit.fit) < 21
                                    else fit.fit[self.popsize // 2])
        if len(fit.histbest) > 2e4: # 10 + 30*N/sp.popsize:
            fit.histbest.pop()
            fit.histmedian.pop()
        if len(fit.hist) > 10 + 30*N/sp.popsize:
            fit.hist.pop()


        # --------------------- this one is not used
        if self.opts['CMA_AII']:
            self.aii.tell(solutions, function_values)
            self.flgtelldone = True
            # for output:
            self.mean = self.aii.mean
            self.dC = self.aii.sigmai**2
            self.sigma = self.aii.sigma
            self.D = 1e-11 + (self.aii.r**2)**0.5
            self.more_to_write += [self.aii.sigma_r]
            return

        # --------------------- LOADING THE SOLUTIONS GENOTYPES INTO pop array
        # TODO: clean up inconsistency when an unrepaired solution is available and used
        pop = []  # create pop from input argument solutions || contains the genotypes of the solutions
        for s in solutions:  # use phenotype before Solution.repair()
            if use_sent_solutions: # this is a boolean that seems to be always True
                x = self.sent_solutions.pop(s, None)  # 12.7s vs 11.3s with N,lambda=20,200 || x contains the hash geno, pheno, iterationno
                if x is not None:
                    pop.append(x['geno'])
                    # TODO: keep additional infos or don't pop s from sent_solutions in the first place
                else:
                    # print 'WARNING: solution not found in ``self.sent_solutions`` (is expected for injected solutions)'
                    pop.append(self.gp.geno(s, copy=copy))  # cannot recover the original genotype with boundary handling
                    if check_points in (None, True, 1):
                        self.repair_genotype(pop[-1])  # necessary if pop[-1] was changed or injected by the user.
#            else:  # TODO: to be removed?
#                # print 'WARNING: ``geno`` mapping depreciated'
#                pop.append(self.gp.geno(s, copy=copy))
#                if check_points in (None, True, 1):
#                    self.repair_genotype(pop[-1])  # necessary or not?
#                # print 'repaired'

        mold = self.mean # simply the mean -> [0.1 0.1 0.1 0.1 0.1] for 5D for example
        sigma_fac = 1

        # --------------------- not going in this so far
        # check and normalize each x - m
        # check_points is a flag (None is default: check non-known solutions) or an index list
        # should also a number possible (first check_points points)?
        if check_points not in (None, False, 0, [], ()):  # useful in case of injected solutions and/or adaptive encoding, however is automatic with use_sent_solutions
            try:
                if len(check_points):
                    idx = check_points
            except:
                idx = xrange(sp.popsize)

            for k in idx:
                self.repair_genotype(pop[k])

        # --------------------- SORTING THE REST OF THE SOLUTIONS BY THEIR FITNESS VALUES (the fit.idx found before)
        if type(pop) is not array: # only arrays can be multiple indexed
            pop = array(pop, copy=False)

#        pop = pop[fit.idx]

        orderByConstraint = indeces[1] + indeces[0]
        pop = pop[orderByConstraint]



        # --------------------- the CMA_elitist case, not of concern for now
        if self.opts['CMA_elitist'] and self.best.f < fit.fit[0]:
            if self.best.x_geno is not None:
                xp = [self.best.x_geno]
                # xp = [self.best.xdict['geno']]
                # xp = [self.gp.geno(self.best.x[:])]  # TODO: remove
                # print self.mahalanobisNorm(xp[0]-self.mean)
                self.clip_or_fit_solutions(xp, [0])
                pop = array([xp[0]] + list(pop))
            else:
                print('genotype for elitist not found')

        # --------------------- COMPUTING THE NEW MEAN
        # because the old mean was calculated based on the sp.mu which is fixed we need to remove it
        # self.sp.mean is alway 1.0, is it a scaling factor of some sort?

        #self.mean = mold + self.sp.cmean * (sum(sp.weights * pop[0:sp.mu].T, 1) - mold)

        if sp.wholePopFlag:
            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop.T, 1) - mold)
        else:
            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop[0:sp.mu].T, 1) - mold)
#            self.mean = self.sp.cmean * sum(sp.weights * pop[0:sp.mu].T, 1) # this mean is modified with the +- mold
#            self.mean = mold + self.sp.cmean * (sum(sp.weights * pop.T, 1) - mold)
#        print(">> new mean: "+str(self.mean))


        # check Delta m (this is not default, but could become at some point)
        # CAVE: upper_length=sqrt(2)+2 is too restrictive, test upper_length = sqrt(2*N) thoroughly.
        # simple test case injecting self.mean:
        # self.mean = 1e-4 * self.sigma * np.random.randn(N)
        if 11 < 3 and self.opts['vv'] and check_points:  # TODO: check_points might be an index-list
            cmean = self.sp.cmean / min(1, (sqrt(self.opts['vv']*N)+2) / ( # abuse of cmean
                (sqrt(self.sp.mueff) / self.sp.cmean) *
                self.mahalanobisNorm(self.mean - mold)))
        else:
            cmean = self.sp.cmean

        if 11 < 3:  # plot length of mean - mold
            self.more_to_write += [sqrt(sp.mueff) *
                sum(((1./self.D) * dot(self.B.T, self.mean - mold))**2)**0.5 /
                       self.sigma / sqrt(N) / cmean]

        # get learning rate constants
        cc, c1, cmu = sp.cc, sp.c1, sp.cmu
        if flgseparable:
            cc, c1, cmu = sp.cc_sep, sp.c1_sep, sp.cmu_sep

        # now the real work can start

        # evolution paths
        self.ps = (1-sp.cs) * self.ps + \
                  (sqrt(sp.cs*(2-sp.cs)*sp.mueff)  / self.sigma / cmean) * \
                  dot(self.B, (1./self.D) * dot(self.B.T, (self.mean - mold) / self.sigma_vec))

        # "hsig", correction with self.countiter seems not necessary, also pc starts with zero
        hsig = sum(self.ps**2) / (1-(1-sp.cs)**(2*self.countiter)) / self.N < 2 + 4./(N+1)
        if 11 < 3:
            # hsig = 1
            # sp.cc = 4 / (N + 4)
            # sp.cs = 4 / (N + 4)
            # sp.cc = 1
            # sp.damps = 2  #
            # sp.CMA_on = False
            # c1 = 0  # 2 / ((N + 1.3)**2 + 0 * sp.mu) # 1 / N**2
            # cmu = min([1 - c1, cmu])
            if self.countiter == 1:
                print('parameters modified')
        # hsig = sum(self.ps**2) / self.N < 2 + 4./(N+1)
        # adjust missing variance due to hsig, in 4-D with damps=1e99 and sig0 small
        #       hsig leads to premature convergence of C otherwise
        #hsiga = (1-hsig**2) * c1 * cc * (2-cc)  # to be removed in future
        c1a = c1 - (1-hsig**2) * c1 * cc * (2-cc)  # adjust for variance loss

        if 11 < 3:  # diagnostic data
            self.out['hsigcount'] += 1 - hsig
            if not hsig:
                self.hsiglist.append(self.countiter)
        if 11 < 3:  # diagnostic message
            if not hsig:
                print(str(self.countiter) + ': hsig-stall')
        if 11 < 3:  # for testing purpose
            hsig = 1 # TODO:
            #       put correction term, but how?
            if self.countiter == 1:
                print('hsig=1')

        self.pc = (1-cc) * self.pc + \
                  hsig * (sqrt(cc*(2-cc)*sp.mueff) / self.sigma / cmean) * \
                  (self.mean - mold)  / self.sigma_vec

        # covariance matrix adaptation/udpate
        if sp.CMA_on:
            # assert sp.c1 + sp.cmu < sp.mueff / N  # ??
            assert c1 + cmu <= 1

            # default full matrix case
            if not flgseparable:
                if sp.wholePopFlag:
                    Z = (pop - mold) / (self.sigma * self.sigma_vec)
                else:
                    Z = (pop[0:sp.mu] - mold) / (self.sigma * self.sigma_vec)
                Z = dot((cmu * sp.weights) * Z.T, Z)  # learning rate integrated
                if self.sp.neg.cmuexp:
                    tmp = (pop[-sp.neg.mu:] - mold) / (self.sigma * self.sigma_vec)
                    self.Zneg *= 1 - self.sp.neg.cmuexp  # for some reason necessary?
                    self.Zneg += dot(sp.neg.weights * tmp.T, tmp) - self.C
                    # self.update_exponential(dot(sp.neg.weights * tmp.T, tmp) - 1 * self.C, -1*self.sp.neg.cmuexp)

#                if 11 < 3: # ?3 to 5 times slower??
#                    Z = np.zeros((N,N))
#                    for k in xrange(sp.mu):
#                        z = (pop[k]-mold)
#                        Z += np.outer((cmu * sp.weights[k] / (self.sigma * self.sigma_vec)**2) * z, z)

                self.C *= 1 - c1a - cmu
                self.C += np.outer(c1 * self.pc, self.pc) + Z
#                self.C = np.eye(self.N, self.N);    # NEW , remove it later !!!!!!!! ILYA
                self.dC = np.diag(self.C)  # for output and termination checking

            else: # separable/diagonal linear case
                assert(c1+cmu <= 1)
                Z = np.zeros(N)
                for k in xrange(sp.mu):
                    z = (pop[k]-mold) / (self.sigma * self.sigma_vec) # TODO see above
                    Z += sp.weights[k] * z * z  # is 1-D
                self.C = (1-c1a-cmu) * self.C + c1 * self.pc * self.pc + cmu * Z
                # TODO: self.C *= exp(cmuneg * (N - dot(sp.neg.weights,  **2)
#                self.C = np.eye(self.N, self.N);   # NEW , remove it later !!!!!!!! ILYA
                self.dC = self.C
                self.D = sqrt(self.C)  # C is a 1-D array
                self.itereigenupdated = self.countiter

                # idx = self.mirror_idx_cov()  # take half of mirrored vectors for negative update

        # qqqqqqqqqqq
        if 1 < 3 and np.isfinite(sp.dampsvec):
            if self.countiter == 1:
                print("WARNING: CMA_dampsvec option is experimental")
            sp.dampsvec *= np.exp(sp.dampsvec_fading/self.N)
            # TODO: rank-lambda update: *= (1 + sum(z[z>1]**2-1) * exp(sum(z[z<1]**2-1))
            self.sigma_vec *= np.exp((sp.cs/sp.dampsvec/2) * (self.ps**2 - 1))
            # self.sigma_vec *= np.exp((sp.cs/sp.dampsvec) * (abs(self.ps) - (2/np.pi)**0.5))
            self.more_to_write += [exp(np.mean((self.ps**2 - 1)**2))]
            # TODO: rank-mu update

        # step-size adaptation, adapt sigma
        if 1 < 3:  #
            self.sigma *= sigma_fac * \
                            np.exp((min((1, (sp.cs/sp.damps) *
                                    (sqrt(sum(self.ps**2))/self.const.chiN - 1)))))
        else:
            self.sigma *= sigma_fac * \
                            np.exp((min((1000, (sp.cs/sp.damps/2) *
                                    (sum(self.ps**2)/N - 1)))))
#        if 11 < 3:
#            # derandomized MSR = natural gradient descent using mean(z**2) instead of mu*mean(z)**2
#            lengths = array([sum(z**2)**0.5 for z in self.arz[fit.idx[:self.sp.mu]]])
#            # print lengths[0::int(self.sp.mu/5)]
#            self.sigma *= np.exp(self.sp.mueff**0.5 * dot(self.sp.weights, lengths / self.const.chiN - 1))**(2/(N+1))
#
#        if 11 < 3 and self.opts['vv']:
#            if self.countiter < 2:
#                print('constant sigma applied')
#                print(self.opts['vv'])  # N=10,lam=10: 0.8 is optimal
#            self.sigma = self.opts['vv'] * self.sp.mueff * sum(self.mean**2)**0.5 / N

        if self.sigma * min(self.dC)**0.5 < self.opts['minstd']:
            self.sigma = self.opts['minstd'] / min(self.dC)**0.5
        # g = self.countiter
        # N = self.N
        mindx = eval(self.opts['mindx']) if type(self.opts['mindx']) == type('') else self.opts['mindx']
        if self.sigma * min(self.D) < mindx:  # TODO: sigma_vec is missing here
            self.sigma = mindx / min(self.D)

        if self.sigma > 1e9 * self.sigma0:
            alpha = self.sigma / max(self.D)
            self.multiplyC(alpha)
            self.sigma /= alpha**0.5
            self.opts['tolupsigma'] /= alpha**0.5  # to be compared with sigma

        # TODO increase sigma in case of a plateau?

        # Uncertainty noise measurement is done on an upper level

#        # output, has moved up, e.g. as part of fmin, TODO to be removed
#        if 11 < 3 and self.opts['verb_log'] > 0 and (self.countiter < 4 or
#                                          self.countiter % self.opts['verb_log'] == 0):
#            # this assumes that two logger with the same name access the same data!
#            CMADataLogger(self.opts['verb_filenameprefix']).register(self, append=True).add()
#            # self.writeOutput(solutions[fit.idx[0]])

        self.flgtelldone = True

    def result(self):
        """return ``(xbest, f(xbest), evaluations_xbest, evaluations, iterations, pheno(xmean), effective_stds)``"""
        # TODO: how about xcurrent?
        return self.best.get() + (
            self.countevals, self.countiter, self.gp.pheno(self.mean), self.gp.scales * self.sigma * self.sigma_vec * self.dC**0.5)

    def clip_or_fit_solutions(self, pop, idx):
        """make sure that solutions fit to sample distribution, this interface will probably change.

        In particular the frequency of long vectors appearing in pop[idx] - self.mean is limited.

        """
        for k in idx:
            self.repair_genotype(pop[k])

    def repair_genotype(self, x):
        """make sure that solutions fit to sample distribution, this interface will probably change.

        In particular the frequency of x - self.mean being long is limited.

        """
        mold = self.mean
        if 1 < 3:  # hard clip at upper_length
            upper_length = self.N**0.5 + 2 * self.N / (self.N+2)  # should become an Option, but how? e.g. [0, 2, 2]
            fac = self.mahalanobisNorm(x - mold) / upper_length

            if fac > 1:
                x = (x - mold) / fac + mold
                # print self.countiter, k, fac, self.mahalanobisNorm(pop[k] - mold)
                # adapt also sigma: which are the trust-worthy/injected solutions?
            elif 11 < 3:
                return exp(np.tanh(((upper_length*fac)**2/self.N-1)/2) / 2)
        else:
            if 'checktail' not in self.__dict__:  # hasattr(self, 'checktail')
                raise NotImplementedError
                # from check_tail_smooth import CheckTail  # for the time being
                # self.checktail = CheckTail()
                # print('untested feature checktail is on')
            fac = self.checktail.addchin(self.mahalanobisNorm(x - mold))

            if fac < 1:
                x = fac * (x - mold) + mold

        return 1.0  # sigma_fac, not in use


    #____________________________________________________________
    #____________________________________________________________
    #
    def updateBD(self):
        """update internal variables for sampling the distribution with the
        current covariance matrix C. This method is O(N^3), if C is not diagonal.

        """
        # itereigenupdated is always up-to-date in the diagonal case
        # just double check here
        if self.itereigenupdated == self.countiter:
            return

        if self.sp.neg.cmuexp:  # cave:
            self.update_exponential(self.Zneg, -self.sp.neg.cmuexp)
            # self.C += self.Zpos  # pos update after Zneg would be the correct update, overall:
            # self.C = self.Zpos + Cs * Mh.expms(-self.sp.neg.cmuexp*Csi*self.Zneg*Csi) * Cs
            self.Zneg = np.zeros((self.N, self.N))

        if self.sigma_vec is not 1 and not np.all(self.sigma_vec == 1):
            self.C = dot(dot(np.diag(self.sigma_vec), self.C), np.diag(self.sigma_vec))
            self.sigma_vec[:] = 1

        if self.opts['CMA_const_trace'] in (True, 1, 2):  # normalize trace of C
            if self.opts['CMA_const_trace'] == 2:
                s = np.exp(np.mean(np.log(self.dC)))
            else:
                s = np.mean(self.dC)
            self.C /= s
            self.dC /= s
        self.C = (self.C + self.C.T) / 2
        # self.C = np.triu(self.C) + np.triu(self.C,1).T  # should work as well
        # self.D, self.B = eigh(self.C) # hermitian, ie symmetric C is assumed

        if type(self.opts['CMA_eigenmethod']) == type(1):
            print('WARNING: option CMA_eigenmethod should be a function, not an integer')
            if self.opts['CMA_eigenmethod'] == -1:
                # pygsl
                # easy to install (well, in Windows install gsl binaries first,
                # set system path to respective libgsl-0.dll (or cp the dll to
                # python\DLLS ?), in unzipped pygsl edit
                # gsl_dist/gsl_site_example.py into gsl_dist/gsl_site.py
                # and run "python setup.py build" and "python setup.py install"
                # in MINGW32)
                if 1 < 3:  # import pygsl on the fly
                    try:
                        import pygsl.eigen.eigenvectors  # TODO efficient enough?
                    except ImportError:
                        print('WARNING: could not find pygsl.eigen module, either install pygsl \n' +
                              '  or set option CMA_eigenmethod=1 (is much slower), option set to 1')
                        self.opts['CMA_eigenmethod'] = 0  # use 0 if 1 is too slow

                    self.D, self.B = pygsl.eigen.eigenvectors(self.C)

            elif self.opts['CMA_eigenmethod'] == 0:
                # TODO: thoroughly test np.linalg.eigh
                #       numpy.linalg.eig crashes in 200-D
                #       and EVecs with same EVals are not orthogonal
                self.D, self.B = np.linalg.eigh(self.C)  # self.B[i] is a row and not an eigenvector
            else:  # is overall two;ten times slower in 10;20-D
                self.D, self.B = Misc.eig(self.C)  # def eig, see below
        else:
            self.D, self.B = self.opts['CMA_eigenmethod'](self.C)


        # assert(sum(self.D-DD) < 1e-6)
        # assert(sum(sum(np.dot(BB, BB.T)-np.eye(self.N))) < 1e-6)
        # assert(sum(sum(np.dot(BB * DD, BB.T) - self.C)) < 1e-6)
        idx = np.argsort(self.D)
        self.D = self.D[idx]
        self.B = self.B[:,idx]  # self.B[i] is a row, columns self.B[:,i] are eigenvectors
        # assert(all(self.B[self.countiter % self.N] == self.B[self.countiter % self.N,:]))

        # qqqqqqqqqq
        if 11 < 3:  # limit condition number to 1e13
            climit = 1e13  # cave: conditioncov termination is 1e14
            if self.D[-1] / self.D[0] > climit:
                self.D += self.D[-1] / climit
            for i in xrange(self.N):
                self.C[i][i] += self.D[-1] / climit

        if 11 < 3 and any(abs(sum(self.B[:,0:self.N-1] * self.B[:,1:], 0)) > 1e-6):
            print('B is not orthogonal')
            print(self.D)
            print(sum(self.B[:,0:self.N-1] * self.B[:,1:], 0))
        else:
            # is O(N^3)
            # assert(sum(abs(self.C - np.dot(self.D * self.B,  self.B.T))) < N**2*1e-11)
            pass
        self.D **= 0.5
        self.itereigenupdated = self.countiter

    def multiplyC(self, alpha):
        """multiply C with a scalar and update all related internal variables (dC, D,...)"""
        self.C *= alpha
        if self.dC is not self.C:
            self.dC *= alpha
        self.D *= alpha**0.5
    def update_exponential(self, Z, eta, BDpair=None):
        ''''''
        if eta == 0:
            return
        if BDpair:
            B, D = BDpair
        else:
            D, B = self.opts['CMA_eigenmethod'](self.C)
            D **= 0.5
        Csi = dot(B, (B / D).T)
        Cs = dot(B, (B * D).T)
        self.C = dot(Cs, dot(Mh.expms(eta * dot(Csi, dot(Z, Csi)), self.opts['CMA_eigenmethod']), Cs))

    #____________________________________________________________
    #____________________________________________________________
    #
    def _updateCholesky(self, A, Ainv, p, alpha, beta):
        """not yet implemented"""
        # BD is A, p is A*Normal(0,I) distributed
        # input is assumed to be numpy arrays
        # Ainv is needed to compute the evolution path
        # this is a stump and is not tested

        raise _Error("not yet implemented")
        # prepare
        alpha = float(alpha)
        beta = float(beta)
        y = np.dot(Ainv, p)
        y_sum = sum(y**2)

        # compute scalars
        tmp = sqrt(1 + beta * y_sum / alpha)
        fac = (sqrt(alpha) / sum(y**2)) * (tmp - 1)
        facinv = (1. / (sqrt(alpha) * sum(y**2))) * (1 - 1. / tmp)

        # update matrices
        A *= sqrt(alpha)
        A += np.outer(fac * p, y)
        Ainv /= sqrt(alpha)
        Ainv -= np.outer(facinv * y, np.dot(y.T, Ainv))

    #____________________________________________________________
    #____________________________________________________________
    def feedForResume(self, X, function_values):
        ''''''
        if self.countiter > 0:
            print('WARNING: feed should generally be used with a new object instance')
        if len(X) != len(function_values):
            raise _Error('number of solutions ' + str(len(X)) +
                ' and number function values ' +
                str(len(function_values))+' must not differ')
        popsize = self.sp.popsize
        if (len(X) % popsize) != 0:
            raise _Error('number of solutions ' + str(len(X)) +
                    ' must be a multiple of popsize (lambda) ' +
                    str(popsize))
        for i in xrange(len(X) / popsize):
            # feed in chunks of size popsize
            self.ask()  # a fake ask, mainly for a conditioned calling of updateBD
                        # and secondary to get possibly the same random state
            self.tell(X[i*popsize:(i+1)*popsize], function_values[i*popsize:(i+1)*popsize])

    #____________________________________________________________
    #____________________________________________________________
    def readProperties(self):
        """reads dynamic parameters from property file (not implemented)
        """
        print('not yet implemented')

    #____________________________________________________________
    #____________________________________________________________
    def mahalanobisNorm(self, dx):
        ''''''
        return sqrt(sum((self.D**-1 * np.dot(self.B.T, dx))**2)) / self.sigma

    #____________________________________________________________
    #____________________________________________________________
    #
    def timesCroot(self, mat):
        """return C**0.5 times mat, where mat can be a vector or matrix.
        Not functional, because _Croot=C**0.5 is never computed (should be in updateBD)
        """
        print("WARNING: timesCroot is not yet tested")
        if self.opts['CMA_diagonal'] is True \
                       or self.countiter <= self.opts['CMA_diagonal']:
            res = (self._Croot * mat.T).T
        else:
            res = np.dot(self._Croot, mat)
        return res
    def divCroot(self, mat):
        """return C**-1/2 times mat, where mat can be a vector or matrix"""
        print("WARNING: divCroot is not yet tested")
        if self.opts['CMA_diagonal'] is True \
                       or self.countiter <= self.opts['CMA_diagonal']:
            res = (self._Crootinv * mat.T).T
        else:
            res = np.dot(self._Crootinv, mat)
        return res

    #____________________________________________________________
    #____________________________________________________________
    def disp_annotation(self):
        """print annotation for `disp()`"""
        print('Iterat #Fevals   function value     axis ratio  sigma   minstd maxstd min:sec')
        sys.stdout.flush()

    #____________________________________________________________
    #____________________________________________________________
    def disp(self, modulo=None):  # TODO: rather assign opt['verb_disp'] as default?
        """prints some infos according to `disp_annotation()`, if
        ``iteration_counter % modulo == 0``

        """
        if modulo is None:
            modulo = self.opts['verb_disp']

        # console display
        if modulo:
            if (self.countiter-1) % (10 * modulo) < 1:
                self.disp_annotation()
            if self.countiter > 0 and (self.stop() or self.countiter < 4
                              or self.countiter % modulo < 1):
                if self.opts['verb_time']:
                    toc = self.elapsed_time()
                    stime = str(int(toc//60))+':'+str(round(toc%60,1))
                else:
                    stime = ''
                # __ANNOT__
                print(' '.join((repr(self.countiter).rjust(5),
                                repr(self.countevals+int(self.unsuccess)).rjust(7),
                                '%.15e' % (min(self.fit.fit)),
                                '%4.1e' % (self.D.max()/self.D.min()),
                                '%6.2e' % self.sigma,
                                '%6.0e' % (self.sigma * sqrt(min(self.dC))),
                                '%6.0e' % (self.sigma * sqrt(max(self.dC))),
                                stime)))
                # if self.countiter < 4:
                sys.stdout.flush()

class Options(dict):

    # @classmethod # self is the class, not the instance
    # @property
    # def default(self):
    #     """returns all options with defaults"""
    #     return fmin([],[])

    @staticmethod
    def defaults():
        """return a dictionary with default option values and description,
        calls `fmin([], [])`"""
        return fmin([], [])

    @staticmethod
    def versatileOptions():
        """return list of options that can be changed at any time (not only be
        initialized), however the list might not be entirely up to date. The
        string ' #v ' in the default value indicates a 'versatile' option
        that can be changed any time.

        """
        return tuple(sorted(i[0] for i in list(Options.defaults().items()) if i[1].find(' #v ') > 0))

    def __init__(self, s=None, unchecked=False):
        """return an `Options` instance, either with the default options,
        if ``s is None``, or with all options whose name or description
        contains `s`, if `s` is a string (case is disregarded),
        or with entries from dictionary `s` as options, not complemented
        with default options or settings

        Returns: see above.

        """
        # if not Options.defaults:  # this is different from self.defaults!!!
        #     Options.defaults = fmin([],[])
        if s is None:
            super(Options, self).__init__(Options.defaults())
            # self = Options.defaults()
        elif type(s) is str:
            super(Options, self).__init__(Options().match(s))
            # we could return here
        else:
            super(Options, self).__init__(s)

        if not unchecked:
            for key in list(self.keys()):
                if key not in Options.defaults():
                    print('Warning in cma.Options.__init__(): invalid key ``' + str(key) + '`` popped')
                    self.pop(key)
        # self.evaluated = False  # would become an option entry

    def init(self, dict_or_str, val=None, warn=True):

        #dic = dict_or_key if val is None else {dict_or_key:val}
        dic = dict_or_str
        if val is not None:
            dic = {dict_or_str:val}

        for key, val in list(dic.items()):
            if key not in Options.defaults():
                # TODO: find a better solution?
                if warn:
                    print('Warning in cma.Options.init(): key ' +
                        str(key) + ' ignored')
            else:
                self[key] = val

        return self

    def set(self, dic, val=None, warn=True):

        if val is not None:  # dic is a key in this case
            dic = {dic:val}  # compose a dictionary
        for key, val in list(dic.items()):
            if key in Options.versatileOptions():
                self[key] = val
            elif warn:
                print('Warning in cma.Options.set(): key ' + str(key) + ' ignored')
        return self  # to allow o = Options(o).set(new)

    def complement(self):
        """add all missing options with their default values"""

        for key in Options.defaults():
            if key not in self:
                self[key] = Options.defaults()[key]
        return self

    def settable(self):
        """return the subset of those options that are settable at any
        time.

        Settable options are in `versatileOptions()`, but the
        list might be incomlete.

        """
        return Options([i for i in list(self.items())
                                if i[0] in Options.versatileOptions()])

    def __call__(self, key, default=None, loc=None):

        try:
            val = self[key]
        except:
            return self.match(key)

        if loc is None:
            loc = self  # TODO: this hack is not so useful: popsize could be there, but N is missing
        try:
            if type(val) is str:
                val = val.split('#')[0].strip()  # remove comments
                if type(val) == type('') and key.find('filename') < 0 and key.find('mindx') < 0:
                    val = eval(val, globals(), loc)
            # invoke default
            # TODO: val in ... fails with array type, because it is applied element wise!
            # elif val in (None,(),[],{}) and default is not None:
            elif val is None and default is not None:
                val = eval(str(default), globals(), loc)
        except:
            pass  # slighly optimistic: the previous is bug-free
        return val

    def eval(self, key, default=None, loc=None):

        self[key] = self(key, default, loc)
        return self[key]

    def evalall(self, loc=None):
        """Evaluates all option values in environment `loc`.

        :See: `eval()`

        """
        # TODO: this needs rather the parameter N instead of loc
        if 'N' in list(loc.keys()):  # TODO: __init__ of CMA can be simplified
            popsize = self('popsize', Options.defaults()['popsize'], loc)
            for k in list(self.keys()):
                self.eval(k, Options.defaults()[k],
                          {'N':loc['N'], 'popsize':popsize})
        return self

    def match(self, s=''):
        """return all options that match, in the name or the description,
        with string `s`, case is disregarded.

        Example: ``cma.Options().match('verb')`` returns the verbosity options.

        """
        match = s.lower()
        res = {}
        for k in sorted(self):
            s = str(k) + '=\'' + str(self[k]) + '\''
            if match in s.lower():
                res[k] = self[k]
        return Options(res)

    def pp(self):
        pprint(self)

    def printme(self, linebreak=80):
        for i in sorted(Options.defaults().items()):
            s = str(i[0]) + "='" + str(i[1]) + "'"
            a = s.split(' ')

            # print s in chunks
            l = ''  # start entire to the left
            while a:
                while a and len(l) + len(a[0]) < linebreak:
                    l += ' ' + a.pop(0)
                print(l)
                l = '        '  # tab for subsequent lines

class CMAStopDict(dict):
    """keep and update a termination condition dictionary, which is
    "usually" empty and returned by `CMAEvolutionStrategy.stop()`.

    Details
    -------
    This could be a nested class, but nested classes cannot be serialized.

    :See: `stop()`

    """
    def __init__(self, d={}):
        update = (type(d) == CMAEvolutionStrategy)
        inherit = (type(d) == CMAStopDict)
        super(CMAStopDict, self).__init__({} if update else d)
        self._stoplist = d._stoplist if inherit else []    # multiple entries
        self.lastiter = d.lastiter if inherit else 0  # probably not necessary
        if update:
            self._update(d)
            
        

    def __call__(self, es):
        """update the dictionary"""
        return self._update(es)

    def _addstop(self, key, cond, val=None):
        if cond:
            self.stoplist.append(key)  # can have the same key twice
            if key in list(self.opts.keys()):
                val = self.opts[key]
            self[key] = val

    def _update(self, es):
        """Test termination criteria and update dictionary.

        """
        if es.countiter == self.lastiter:
            if es.countiter == 0:
                self.__init__()
                return self
            try:
                if es == self.es:
                    return self
            except: # self.es not yet assigned
                pass

        self.lastiter = es.countiter
        self.es = es

        self.stoplist = []

        N = es.N
        opts = es.opts
        self.opts = opts  # a hack to get _addstop going

        # fitness: generic criterion, user defined w/o default
        self._addstop('ftarget',
                     es.best.f < opts['ftarget'])
        # maxiter, maxfevals: generic criteria
        self._addstop('maxfevals',
                     es.countevals - 1 >= opts['maxfevals'])
        self._addstop('maxiter',
                     es.countiter >= opts['maxiter'])
        # tolx, tolfacupx: generic criteria
        # tolfun, tolfunhist (CEC:tolfun includes hist)
        
        
        self._addstop('tolx',
                     all([es.sigma*xi < opts['tolx'] for xi in es.pc]) and \
                     all([es.sigma*xi < opts['tolx'] for xi in sqrt(es.dC)]))
        

        self._addstop('tolfacupx',
                     any([es.sigma * sig > es.sigma0 * opts['tolfacupx']
                          for sig in sqrt(es.dC)]))
        self._addstop('tolfun',
                     es.fit.fit[-1] - es.fit.fit[0] < opts['tolfun'] and \
                     max(es.fit.hist) - min(es.fit.hist) < opts['tolfun'])
        self._addstop('tolfunhist',
                     len(es.fit.hist) > 9 and \
                     max(es.fit.hist) - min(es.fit.hist) <  opts['tolfunhist'])

        # worst seen false positive: table N=80,lam=80, getting worse for fevals=35e3 \approx 50 * N**1.5
        # but the median is not so much getting worse
        # / 5 reflects the sparsity of histbest/median
        # / 2 reflects the left and right part to be compared
        l = int(max(opts['tolstagnation'] / 5. / 2, len(es.fit.histbest) / 10));
        # TODO: why max(..., len(histbest)/10) ???
        # TODO: the problem in the beginning is only with best ==> ???
        if 11 < 3:  #
            print(es.countiter, (opts['tolstagnation'], es.countiter > N * (5 + 100 / es.popsize),
                        len(es.fit.histbest) > 100,
                        np.median(es.fit.histmedian[:l]) >= np.median(es.fit.histmedian[l:2*l]),
                        np.median(es.fit.histbest[:l]) >= np.median(es.fit.histbest[l:2*l])))
        # equality should handle flat fitness
        self._addstop('tolstagnation', # leads sometimes early stop on ftablet, fcigtab, N>=50?
                    1 < 3 and opts['tolstagnation'] and es.countiter > N * (5 + 100 / es.popsize) and
                    len(es.fit.histbest) > 100 and 2*l < len(es.fit.histbest) and
                    np.median(es.fit.histmedian[:l]) >= np.median(es.fit.histmedian[l:2*l]) and
                    np.median(es.fit.histbest[:l]) >= np.median(es.fit.histbest[l:2*l]))
        # iiinteger: stagnation termination can prevent to find the optimum

        self._addstop('tolupsigma', opts['tolupsigma'] and
                      es.sigma / es.sigma0 / np.max(es.D) > opts['tolupsigma'])

        if 11 < 3 and 2*l < len(es.fit.histbest):  # TODO: this might go wrong, because the nb of written columns changes
            tmp = np.array((-np.median(es.fit.histmedian[:l]) + np.median(es.fit.histmedian[l:2*l]),
                        -np.median(es.fit.histbest[:l]) + np.median(es.fit.histbest[l:2*l])))
            es.more_to_write += [(10**t if t < 0 else t + 1) for t in tmp] # the latter to get monotonicy

        if 1 < 3:
            # non-user defined, method specific
            # noeffectaxis (CEC: 0.1sigma), noeffectcoord (CEC:0.2sigma), conditioncov
            self._addstop('noeffectcoord',
                         any([es.mean[i] == es.mean[i] + 0.2*es.sigma*sqrt(es.dC[i])
                              for i in xrange(N)]))
            if opts['CMA_diagonal'] is not True and es.countiter > opts['CMA_diagonal']:
                i = es.countiter % N
                self._addstop('noeffectaxis',
                             sum(es.mean == es.mean + 0.1 * es.sigma * es.D[i] * es.B[:, i]) == N)
            self._addstop('conditioncov',
                         es.D[-1] > 1e7 * es.D[0], 1e14)  # TODO

            self._addstop('callback', es.callbackstop)  # termination_callback
        if len(self):
            self._addstop('flat fitness: please (re)consider how to compute the fitness more elaborate',
                         len(es.fit.hist) > 9 and \
                         max(es.fit.hist) == min(es.fit.hist))
        if 11 < 3 and opts['vv'] == 321:
            self._addstop('||xmean||^2<ftarget', sum(es.mean**2) <= opts['ftarget'])

        return self

class BaseDataLogger(object):
    """"abstract" base class for a data logger that can be used with an `OOOptimizer`"""
    def add(self, optim=None, more_data=[]):
        """abstract method, add a "data point" from the state of `optim` into the
        logger, the argument `optim` can be omitted if it was `register()`-ed before,
        acts like an event handler"""
        raise NotImplementedError()
    def register(self, optim):
        """abstract method, register an optimizer `optim`, only needed if `add()` is
        called without a value for the `optim` argument"""
        self.optim = optim
    def disp(self):
        """display some data trace (not implemented)"""
        print('method BaseDataLogger.disp() not implemented, to be done in subclass ' + str(type(self)))
    def plot(self):
        """plot data (not implemented)"""
        print('method BaseDataLogger.plot() is not implemented, to be done in subclass ' + str(type(self)))
    def data(self):
        """return logged data in a dictionary (not implemented)"""
        print('method BaseDataLogger.data() is not implemented, to be done in subclass ' + str(type(self)))

class CMADataLogger(BaseDataLogger):  # might become a dict at some point

    default_prefix = 'outcmaes'
    # names = ('axlen','fit','stddev','xmean','xrecentbest')
    # key_names_with_annotation = ('std', 'xmean', 'xrecent')

    def __init__(self, name_prefix=default_prefix, modulo=1, append=False):
        """initialize logging of data from a `CMAEvolutionStrategy` instance,
        default modulo expands to 1 == log with each call

        """
        # super(CMAData, self).__init__({'iter':[], 'stds':[], 'D':[], 'sig':[], 'fit':[], 'xm':[]})
        # class properties:
        self.file_names = ('axlen','fit','stddev','xmean','xrecentbest') # used in load, however hard-coded in add
        self.key_names = ('D', 'f', 'std', 'xmean', 'xrecent') # used in load, however hard-coded in plot
        self.key_names_with_annotation = ('std', 'xmean', 'xrecent') # used in load
        self.modulo = modulo  # allows calling with None
        self.append = append
        self.counter = 0  # number of calls of add, should initial value depend on `append`?
        self.name_prefix = name_prefix if name_prefix else CMADataLogger.default_prefix
        if type(self.name_prefix) == CMAEvolutionStrategy:
            self.name_prefix = self.name_prefix.opts.eval('verb_filenameprefix')
        self.registered = False

    def register(self, es, append=None, modulo=None):
        """register a `CMAEvolutionStrategy` instance for logging,
        ``append=True`` appends to previous data logged under the same name,
        by default previous data are overwritten.
        """
        if type(es) != CMAEvolutionStrategy:
            raise TypeError("only class CMAEvolutionStrategy can be registered for logging")
        self.es = es
        if append is not None:
            self.append = append
        if modulo is not None:
            self.modulo = modulo
        if not self.append and self.modulo != 0:
            self.initialize()  # write file headers
        self.registered = True
        return self

    def initialize(self, modulo=None):
        """reset logger, overwrite original files, `modulo`: log only every modulo call"""
        if modulo is not None:
            self.modulo = modulo
        try:
            es = self.es  # must have been registered
        except AttributeError:
            pass  # TODO: revise usage of es... that this can pass
            raise _Error('call register() before initialize()')

        self.counter = 0  # number of calls of add

        # write headers for output
        fn = self.name_prefix + 'fit.dat'
        strseedtime = 'seed=%d, %s' % (es.opts['seed'], time.asctime())

        try:
            with open(fn, 'w') as f:
                f.write('% # columns="iteration, evaluation, sigma, axis ratio, ' +
                        'bestever, best, median, worst objective function value, ' +
                        'further objective values of best", ' +
                        strseedtime +
                        # strftime("%Y/%m/%d %H:%M:%S", localtime()) + # just asctime() would do
                        '\n')
        except (IOError, OSError):
            print('could not open file ' + fn)

        fn = self.name_prefix + 'axlen.dat'
        try:
            f = open(fn, 'w')
            f.write('%  columns="iteration, evaluation, sigma, max axis length, ' +
                    ' min axis length, all principle axes lengths ' +
                    ' (sorted square roots of eigenvalues of C)", ' +
                    strseedtime +
                    '\n')
            f.close()
        except (IOError, OSError):
            print('could not open file ' + fn)
        finally:
            f.close()
        fn = self.name_prefix + 'stddev.dat'
        try:
            f = open(fn, 'w')
            f.write('% # columns=["iteration, evaluation, sigma, void, void, ' +
                    ' stds==sigma*sqrt(diag(C))", ' +
                    strseedtime +
                    '\n')
            f.close()
        except (IOError, OSError):
            print('could not open file ' + fn)
        finally:
            f.close()

        fn = self.name_prefix + 'xmean.dat'
        try:
            with open(fn, 'w') as f:
                f.write('% # columns="iteration, evaluation, void, void, void, xmean", ' +
                        strseedtime)
                f.write(' # scaling_of_variables: ')
                if np.size(es.gp.scales) > 1:
                    f.write(' '.join(map(str, es.gp.scales)))
                else:
                    f.write(str(es.gp.scales))
                f.write(', typical_x: ')
                if np.size(es.gp.typical_x) > 1:
                    f.write(' '.join(map(str, es.gp.typical_x)))
                else:
                    f.write(str(es.gp.typical_x))
                f.write('\n')
                f.close()
        except (IOError, OSError):
            print('could not open/write file ' + fn)

        fn = self.name_prefix + 'xrecentbest.dat'
        try:
            with open(fn, 'w') as f:
                f.write('% # iter+eval+sigma+0+fitness+xbest, ' +
                        strseedtime +
                        '\n')
        except (IOError, OSError):
            print('could not open/write file ' + fn)

        return self
    # end def __init__

    def load(self, filenameprefix=None):
        """loads data from files written and return a data dictionary, *not*
        a prerequisite for using `plot()` or `disp()`.

        Argument `filenameprefix` is the filename prefix of data to be loaded (five files),
        by default ``'outcmaes'``.

        Return data dictionary with keys `xrecent`, `xmean`, `f`, `D`, `std`

        """
        if not filenameprefix:
            filenameprefix = self.name_prefix
        for i in xrange(len(self.file_names)):
            fn = filenameprefix + self.file_names[i] + '.dat'
            try:
                self.__dict__[self.key_names[i]] = _fileToMatrix(fn)
            except:
                print('WARNING: reading from file "' + fn + '" failed')
            if self.key_names[i] in self.key_names_with_annotation:
                self.__dict__[self.key_names[i]].append(self.__dict__[self.key_names[i]][-1])  # copy last row to later fill in annotation position for display
            self.__dict__[self.key_names[i]] = array(self.__dict__[self.key_names[i]], copy=False)
        return self

    def add(self, es=None, more_data=[], modulo=None): # TODO: find a different way to communicate current x and f
        """append some logging data from `CMAEvolutionStrategy` class instance `es`,
        if ``number_of_times_called % modulo`` equals to zero, never if ``modulo==0``.

        The sequence ``more_data`` must always have the same length.

        When used for a different optimizer class, this function can be
        (easily?) adapted by changing the assignments under INTERFACE
        in the implemention.

        """
        self.counter += 1
        mod = modulo if modulo is not None else self.modulo
        if mod == 0 or (self.counter > 3 and self.counter % mod):
            return
        if es is None:
            try:
                es = self.es  # must have been registered
            except AttributeError :
                raise _Error('call `add` with argument `es` or ``register(es)`` before ``add()``')
        elif not self.registered:
            self.register(es) # calls initialize

        # --- INTERFACE, can be changed if necessary ---
        if type(es) is not CMAEvolutionStrategy: # not necessary
            print('WARNING: <type \'CMAEvolutionStrategy\'> expected, found '
                            + str(type(es)) + ' in method CMADataLogger.add')
        evals = es.countevals
        iteration = es.countiter
        sigma = es.sigma
        axratio = es.D.max()/es.D.min()
        xmean = es.mean # TODO: should be optionally phenotype?
        fmean_noise_free = es.fmean_noise_free
        fmean = es.fmean
        try:
            besteverf = es.best.f
            bestf = es.fit.fit[0]
            medianf = es.fit.fit[es.sp.popsize//2]
            worstf = es.fit.fit[-1]
        except:
            if self.counter > 1: # first call without f-values is OK
                raise
        try:
            xrecent = es.best.last.x
        except:
            xrecent = None
        maxD = es.D.max()
        minD = es.D.min()
        diagD = es.D
        diagC = es.sigma*es.sigma_vec*sqrt(es.dC)
        more_to_write = es.more_to_write
        es.more_to_write = []
        # --- end interface ---

        try:

            # fit
            if self.counter > 1:
                fn = self.name_prefix + 'fit.dat'
                with open(fn, 'a') as f:
                    f.write(str(iteration) + ' '
                            + str(evals) + ' '
                            + str(sigma) + ' '
                            + str(axratio) + ' '
                            + str(besteverf) + ' '
                            + '%.16e' % bestf + ' '
                            + str(medianf) + ' '
                            + str(worstf) + ' '
                            # + str(es.sp.popsize) + ' '
                            # + str(10**es.noiseS) + ' '
                            # + str(es.sp.cmean) + ' '
                            + ' '.join(str(i) for i in more_to_write)
                            + ' '.join(str(i) for i in more_data)
                            + '\n')
            # axlen
            fn = self.name_prefix + 'axlen.dat'
            with open(fn, 'a') as f:  # does not rely on reference counting
                f.write(str(iteration) + ' '
                        + str(evals) + ' '
                        + str(sigma) + ' '
                        + str(maxD) + ' '
                        + str(minD) + ' '
                        + ' '.join(map(str, diagD))
                        + '\n')
            # stddev
            fn = self.name_prefix + 'stddev.dat'
            with open(fn, 'a') as f:
                f.write(str(iteration) + ' '
                        + str(evals) + ' '
                        + str(sigma) + ' '
                        + '0 0 '
                        + ' '.join(map(str, diagC))
                        + '\n')
            # xmean
            fn = self.name_prefix + 'xmean.dat'
            with open(fn, 'a') as f:
                if iteration < 1: # before first iteration
                    f.write('0 0 0 0 0 '
                            + ' '.join(map(str, xmean))
                            + '\n')
                else:
                    f.write(str(iteration) + ' '
                            + str(evals) + ' '
                            # + str(sigma) + ' '
                            + '0 '
                            + str(fmean_noise_free) + ' '
                            + str(fmean) + ' '  # TODO: this does not make sense
                            # TODO should be optional the phenotyp?
                            + ' '.join(map(str, xmean))
                            + '\n')
            # xrecent
            fn = self.name_prefix + 'xrecentbest.dat'
            if iteration > 0 and xrecent is not None:
                with open(fn, 'a') as f:
                    f.write(str(iteration) + ' '
                            + str(evals) + ' '
                            + str(sigma) + ' '
                            + '0 '
                            + str(bestf) + ' '
                            + ' '.join(map(str, xrecent))
                            + '\n')

        except (IOError, OSError):
            if iteration <= 1:
                print('could not open/write file')

    def closefig(self):
        pylab.close(self.fighandle)

    def save(self, nameprefix, switch=False):
        """saves logger data to a different set of files, for
        ``switch=True`` also the loggers name prefix is switched to
        the new value

        """
        if not nameprefix or type(nameprefix) is not str:
            _Error('filename prefix must be a nonempty string')

        if nameprefix == self.default_prefix:
            _Error('cannot save to default name "' + nameprefix + '...", chose another name')

        if nameprefix == self.name_prefix:
            return

        for name in CMADataLogger.names:
            open(nameprefix+name+'.dat', 'w').write(open(self.name_prefix+name+'.dat').read())

        if switch:
            self.name_prefix = nameprefix

    def plot(self, fig=None, iabscissa=1, iteridx=None, plot_mean=True,  # TODO: plot_mean default should be False
             foffset=1e-19, x_opt = None, fontsize=10):

        dat = self.load(self.name_prefix)

        try:
            # pylab: prodedural interface for matplotlib
            from  matplotlib.pylab import figure, ioff, ion, subplot, semilogy, hold, plot, grid, \
                 axis, title, text, xlabel, isinteractive, draw, gcf

        except ImportError:
            ImportError('could not find matplotlib.pylab module, function plot() is not available')
            return

        if fontsize and pylab.rcParams['font.size'] != fontsize:
            print('global variable pylab.rcParams[\'font.size\'] set (from ' +
                  str(pylab.rcParams['font.size']) + ') to ' + str(fontsize))
            pylab.rcParams['font.size'] = fontsize  # subtracted in the end, but return can happen inbetween

        if fig:
            figure(fig)
        else:
            figure(325)
            # show()  # should not be necessary
        self.fighandle = gcf()  # fighandle.number

        if iabscissa not in (0,1):
            iabscissa = 1
        interactive_status = isinteractive()
        ioff() # prevents immediate drawing

        dat.x = dat.xmean    # this is the genotyp
        if not plot_mean:
            try:
                dat.x = dat.xrecent
            except:
                pass
        if len(dat.x) < 2:
            print('not enough data to plot')
            return {}

        if iteridx is not None:
            dat.f = dat.f[np.where([x in iteridx for x in dat.f[:,0]])[0],:]
            dat.D = dat.D[np.where([x in iteridx for x in dat.D[:,0]])[0],:]
            iteridx.append(dat.x[-1,1])  # last entry is artificial
            dat.x = dat.x[np.where([x in iteridx for x in dat.x[:,0]])[0],:]
            dat.std = dat.std[np.where([x in iteridx for x in dat.std[:,0]])[0],:]

        if iabscissa == 0:
            xlab = 'iterations'
        elif iabscissa == 1:
            xlab = 'function evaluations'

        # use fake last entry in x and std for line extension-annotation
        if dat.x.shape[1] < 100:
            minxend = int(1.06*dat.x[-2, iabscissa])
            # write y-values for individual annotation into dat.x
            dat.x[-1, iabscissa] = minxend  # TODO: should be ax[1]
            idx = np.argsort(dat.x[-2,5:])
            idx2 = np.argsort(idx)
            if x_opt is None:
                dat.x[-1,5+idx] = np.linspace(np.min(dat.x[:,5:]),
                            np.max(dat.x[:,5:]), dat.x.shape[1]-5)
            else:
                dat.x[-1,5+idx] = np.logspace(np.log10(np.min(abs(dat.x[:,5:]))),
                            np.log10(np.max(abs(dat.x[:,5:]))), dat.x.shape[1]-5)
        else:
            minxend = 0

        if len(dat.f) == 0:
            print('nothing to plot')
            return

        # not in use anymore, see formatter above
        # xticklocs = np.arange(5) * np.round(minxend/4., -int(np.log10(minxend/4.)))

        # dfit(dfit<1e-98) = NaN;

        ioff() # turns update off

        # TODO: if abscissa==0 plot in chunks, ie loop over subsets where dat.f[:,0]==countiter is monotonous

        subplot(2,2,1)
        self.plotdivers(dat, iabscissa, foffset)

        # TODO: modularize also the remaining subplots
        subplot(2,2,2)
        hold(False)
        if x_opt is not None:  # TODO: differentate neg and pos?
            semilogy(dat.x[:, iabscissa], abs(dat.x[:,5:]) - x_opt, '-')
        else:
            plot(dat.x[:, iabscissa], dat.x[:,5:],'-')
        hold(True)
        grid(True)
        ax = array(axis())
        # ax[1] = max(minxend, ax[1])
        axis(ax)
        ax[1] -= 1e-6
        if dat.x.shape[1] < 100:
            yy = np.linspace(ax[2]+1e-6, ax[3]-1e-6, dat.x.shape[1]-5)
            #yyl = np.sort(dat.x[-1,5:])
            idx = np.argsort(dat.x[-1,5:])
            idx2 = np.argsort(idx)
            if x_opt is not None:
                semilogy([dat.x[-1, iabscissa], ax[1]], [abs(dat.x[-1,5:]), yy[idx2]], 'k-') # line from last data point
                semilogy(np.dot(dat.x[-2, iabscissa],[1,1]), array([ax[2]+1e-6, ax[3]-1e-6]), 'k-')
            else:
                # plot([dat.x[-1, iabscissa], ax[1]], [dat.x[-1,5:], yy[idx2]], 'k-') # line from last data point
                plot(np.dot(dat.x[-2, iabscissa],[1,1]), array([ax[2]+1e-6, ax[3]-1e-6]), 'k-')
            # plot(array([dat.x[-1, iabscissa], ax[1]]),
            #      reshape(array([dat.x[-1,5:], yy[idx2]]).flatten(), (2,4)), '-k')
            for i in range(len(idx)):
                # TODOqqq: annotate phenotypic value!?
                # text(ax[1], yy[i], 'x(' + str(idx[i]) + ')=' + str(dat.x[-2,5+idx[i]]))
                text(dat.x[-1,iabscissa], dat.x[-1,5+i], 'x(' + str(i) + ')=' + str(dat.x[-2,5+i]))

        i = 2  # find smallest i where iteration count differs (in case the same row appears twice)
        while i < len(dat.f) and dat.f[-i][0] == dat.f[-1][0]:
            i += 1
        title('Object Variables (' + ('mean' if plot_mean else 'curr best') +
                ', ' + str(dat.x.shape[1]-5) + '-D, popsize~' +
                (str(int((dat.f[-1][1] - dat.f[-i][1]) / (dat.f[-1][0] - dat.f[-i][0])))
                    if len(dat.f.T[0]) > 1 and dat.f[-1][0] > dat.f[-i][0] else 'NA')
                + ')')
        # pylab.xticks(xticklocs)

        # Scaling
        subplot(2,2,3)
        hold(False)
        semilogy(dat.D[:, iabscissa], dat.D[:,5:], '-b')
        hold(True)
        grid(True)
        ax = array(axis())
        # ax[1] = max(minxend, ax[1])
        axis(ax)
        title('Scaling (All Main Axes)')
        # pylab.xticks(xticklocs)
        xlabel(xlab)

        # standard deviations
        subplot(2,2,4)
        hold(False)
        # remove sigma from stds (graphs become much better readible)
        dat.std[:,5:] = np.transpose(dat.std[:,5:].T / dat.std[:,2].T)
        # ax = array(axis())
        # ax[1] = max(minxend, ax[1])
        # axis(ax)
        if 1 < 2 and dat.std.shape[1] < 100:
            # use fake last entry in x and std for line extension-annotation
            minxend = int(1.06*dat.x[-2, iabscissa])
            dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1]
            idx = np.argsort(dat.std[-2,5:])
            idx2 = np.argsort(idx)
            dat.std[-1,5+idx] = np.logspace(np.log10(np.min(dat.std[:,5:])),
                            np.log10(np.max(dat.std[:,5:])), dat.std.shape[1]-5)

            dat.std[-1, iabscissa] = minxend  # TODO: should be ax[1]
            yy = np.logspace(np.log10(ax[2]), np.log10(ax[3]), dat.std.shape[1]-5)
            #yyl = np.sort(dat.std[-1,5:])
            idx = np.argsort(dat.std[-1,5:])
            idx2 = np.argsort(idx)
            # plot(np.dot(dat.std[-2, iabscissa],[1,1]), array([ax[2]+1e-6, ax[3]-1e-6]), 'k-') # vertical separator
            # vertical separator
            plot(np.dot(dat.std[-2, iabscissa],[1,1]), array([np.min(dat.std[-2,5:]), np.max(dat.std[-2,5:])]), 'k-')
            hold(True)
            # plot([dat.std[-1, iabscissa], ax[1]], [dat.std[-1,5:], yy[idx2]], 'k-') # line from last data point
            for i in xrange(len(idx)):
                # text(ax[1], yy[i], ' '+str(idx[i]))
                text(dat.std[-1, iabscissa], dat.std[-1, 5+i], ' '+str(i))
        semilogy(dat.std[:, iabscissa], dat.std[:,5:], '-')
        grid(True)
        title('Standard Deviations in All Coordinates')
        # pylab.xticks(xticklocs)
        xlabel(xlab)
        draw()  # does not suffice
        if interactive_status:
            ion()  # turns interactive mode on (again)
            draw()
        show()

        return self


    #____________________________________________________________
    #____________________________________________________________
    #
    @staticmethod
    def plotdivers(dat, iabscissa, foffset):

        from  matplotlib.pylab import semilogy, hold, grid, \
                 axis, title, text
        fontsize = pylab.rcParams['font.size']

        hold(False)

        dfit = dat.f[:,5]-min(dat.f[:,5])
        dfit[dfit<1e-98] = np.NaN

        if dat.f.shape[1] > 7:
            # semilogy(dat.f[:, iabscissa], abs(dat.f[:,[6, 7, 10, 12]])+foffset,'-k')
            semilogy(dat.f[:, iabscissa], abs(dat.f[:,[6, 7]])+foffset,'-k')
            hold(True)

        # (larger indices): additional fitness data, for example constraints values
        if dat.f.shape[1] > 8:
            # dd = abs(dat.f[:,7:]) + 10*foffset
            # dd = np.where(dat.f[:,7:]==0, np.NaN, dd) # cannot be
            semilogy(dat.f[:, iabscissa], np.abs(dat.f[:,8:]) + 10*foffset, 'm')
            hold(True)

        idx = np.where(dat.f[:,5]>1e-98)[0]  # positive values
        semilogy(dat.f[idx, iabscissa], dat.f[idx,5]+foffset, '.b')
        hold(True)
        grid(True)

        idx = np.where(dat.f[:,5] < -1e-98)  # negative values
        semilogy(dat.f[idx, iabscissa], abs(dat.f[idx,5])+foffset,'.r')

        semilogy(dat.f[:, iabscissa],abs(dat.f[:,5])+foffset,'-b')
        semilogy(dat.f[:, iabscissa], dfit, '-c')

        if 11 < 3:  # delta-fitness as points
            dfit = dat.f[1:, 5] - dat.f[:-1,5]  # should be negative usually
            semilogy(dat.f[1:,iabscissa],  # abs(fit(g) - fit(g-1))
                np.abs(dfit)+foffset, '.c')
            i = dfit > 0
            # print(np.sum(i) / float(len(dat.f[1:,iabscissa])))
            semilogy(dat.f[1:,iabscissa][i],  # abs(fit(g) - fit(g-1))
                np.abs(dfit[i])+foffset, '.r')

        # overall minimum
        i = np.argmin(dat.f[:,5])
        semilogy(dat.f[i, iabscissa]*np.ones(2), dat.f[i,5]*np.ones(2), 'rd')
        # semilogy(dat.f[-1, iabscissa]*np.ones(2), dat.f[-1,4]*np.ones(2), 'rd')

        # AR and sigma
        semilogy(dat.f[:, iabscissa], dat.f[:,3], '-r') # AR
        semilogy(dat.f[:, iabscissa], dat.f[:,2],'-g') # sigma
        semilogy(dat.std[:-1, iabscissa], np.vstack([list(map(max, dat.std[:-1,5:])), list(map(min, dat.std[:-1,5:]))]).T,
                     '-m', linewidth=2)
        text(dat.std[-2, iabscissa], max(dat.std[-2, 5:]), 'max std', fontsize=fontsize)
        text(dat.std[-2, iabscissa], min(dat.std[-2, 5:]), 'min std', fontsize=fontsize)
        ax = array(axis())
        # ax[1] = max(minxend, ax[1])
        axis(ax)
        text(ax[0]+0.01, ax[2], # 10**(log10(ax[2])+0.05*(log10(ax[3])-log10(ax[2]))),
             '.f_recent=' + repr(dat.f[-1,5]) )

        # title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)')
        title('blue:abs(f), cyan:f-min(f), green:sigma, red:axis ratio', fontsize=fontsize-1)
        # pylab.xticks(xticklocs)


    def downsampling(self, factor=10, first=3, switch=True):

        newprefix = self.name_prefix + 'down'
        for name in CMADataLogger.names:
            f = open(newprefix+name+'.dat','w')
            iline = 0
            cwritten = 0
            for line in open(self.name_prefix+name+'.dat'):
                if iline < first or iline % factor == 0:
                    f.write(line)
                    cwritten += 1
                iline += 1
            f.close()
            print('%d' % (cwritten) + ' lines written in ' + newprefix+name+'.dat')
        if switch:
            self.name_prefix += 'down'
        return self

    #____________________________________________________________
    #____________________________________________________________
    #
    def disp(self, idx=100):  # r_[0:5,1e2:1e9:1e2,-10:0]):

        filenameprefix=self.name_prefix

        def printdatarow(dat, iteration):
            """print data of iteration i"""
            i = np.where(dat.f[:, 0] == iteration)[0][0]
            j = np.where(dat.std[:, 0] == iteration)[0][0]
            print('%5d' % (int(dat.f[i,0])) + ' %6d' % (int(dat.f[i,1])) + ' %.14e' % (dat.f[i,5]) +
                  ' %5.1e' % (dat.f[i,3]) +
                  ' %6.2e' % (max(dat.std[j,5:])) + ' %6.2e' % min(dat.std[j,5:]))

        dat = CMADataLogger(filenameprefix).load()
        ndata = dat.f.shape[0]

        # map index to iteration number, is difficult if not all iteration numbers exist
        # idx = idx[np.where(map(lambda x: x in dat.f[:,0], idx))[0]] # TODO: takes pretty long
        # otherwise:
        if idx is None:
            idx = 100
        if np.isscalar(idx):
            # idx = np.arange(0, ndata, idx)
            if idx:
                idx = np.r_[0, 1, idx:ndata-3:idx, -3:0]
            else:
                idx = np.r_[0, 1, -3:0]

        idx = array(idx)
        idx = idx[idx<ndata]
        idx = idx[-idx<=ndata]
        iters = dat.f[idx, 0]
        idxbest = np.argmin(dat.f[:,5])
        iterbest = dat.f[idxbest, 0]

        if len(iters) == 1:
            printdatarow(dat, iters[0])
        else:
            self.disp_header()
            for i in iters:
                printdatarow(dat, i)
            self.disp_header()
            printdatarow(dat, iterbest)
        sys.stdout.flush()
    def disp_header(self):
        heading = 'Iterat Nfevals  function value    axis ratio maxstd   minstd'
        print(heading)

def fmin(func, x0, sigma0=None, args=()
    # the follow string arguments are evaluated, besides the verb_filenameprefix
    , CMA_active='False  # exponential negative update, conducted after the original update'
    , CMA_activefac='1  # learning rate multiplier for active update'
    , CMA_cmean='1  # learning rate for the mean value'
    , CMA_const_trace='False  # normalize trace, value CMA_const_trace=2 normalizes sum log eigenvalues to zero'
    , CMA_diagonal='0*100*N/sqrt(popsize)  # nb of iterations with diagonal covariance matrix, True for always' # TODO 4/ccov_separable?
    , CMA_eigenmethod='np.linalg.eigh  # 0=numpy-s eigh, -1=pygsl, otherwise cma.Misc.eig (slower)'
    , CMA_elitist='False # elitism likely impairs global search performance'
    , CMA_mirrors='popsize < 6  # values <0.5 are interpreted as fraction, values >1 as numbers (rounded), otherwise about 0.16 is used'
    , CMA_mu='None  # parents selection parameter, default is popsize // 2'
    , CMA_on='True  # False or 0 for no adaptation of the covariance matrix'
    , CMA_rankmu='True  # False or 0 for omitting rank-mu update of covariance matrix'
    , CMA_rankmualpha='0.3  # factor of rank-mu update if mu=1, subject to removal, default might change to 0.0'
    , CMA_dampfac='1  #v positive multiplier for step-size damping, 0.3 is close to optimal on the sphere'
    , CMA_dampsvec_fac='np.Inf  # tentative and subject to changes, 0.5 would be a "default" damping for sigma vector update'
    , CMA_dampsvec_fade='0.1  # tentative fading out parameter for sigma vector update'
    , CMA_teststds='None  # factors for non-isotropic initial distr. mainly for test purpose, see scaling_...'
    , CMA_AII='False  # not yet tested'
    , bounds='[None, None]  # lower (=bounds[0]) and upper domain boundaries, each a scalar or a list/vector'
    , eval_parallel='False  # when True, func might be called with more than one solution as first argument'
    , eval_initial_x='False  # '
    , fixed_variables='None  # dictionary with index-value pairs like {0:1.1, 2:0.1} that are not optimized'
    , ftarget='-inf  #v target function value, minimization'
    , incpopsize='2  # in fmin(): multiplier for increasing popsize before each restart'
    , maxfevals='inf  #v maximum number of function evaluations'
    , maxiter='100 + 50 * (N+3)**2 // popsize**0.5  #v maximum number of iterations'
    , mindx='0  #v minimal std in any direction, cave interference with tol*'
    , minstd='0  #v minimal std in any coordinate direction, cave interference with tol*'
    , noise_handling='False  # maximal number of evaluations for noise treatment, only fmin'
    , noise_reevals=' 1.5 + popsize/20  # number of solution to be reevaluated for noise measurement, only fmin'
    , noise_eps='1e-7  # perturbation factor for noise handling reevaluations, only fmin'
    , noise_change_sigma='True  # exponent to default sigma increment'
    , popsize='4+int(3*log(N))  # population size, AKA lambda, number of new solution per iteration'
    , randn='np.random.standard_normal  #v randn((lam, N)) must return an np.array of shape (lam, N)'
    , restarts='0  # in fmin(): number of restarts'
    , restart_from_best='False'
    , scaling_of_variables='None  # scale for each variable, sigma0 is interpreted w.r.t. this scale, in that effective_sigma0 = sigma0*scaling. Internally the variables are divided by scaling_of_variables and sigma is unchanged, default is ones(N)'
    , seed='None  # random number seed'
    , termination_callback='None  #v a function returning True for termination, called after each iteration step and could be abused for side effects'
    , tolfacupx='1e3  #v termination when step-size increases by tolfacupx (diverges). That is, the initial step-size was chosen far too small and better solutions were found far away from the initial solution x0'
    , tolupsigma='1e20  #v sigma/sigma0 > tolupsigma * max(sqrt(eivenvals(C))) indicates "creeping behavior" with usually minor improvements'
    , tolfun='1e-11  #v termination criterion: tolerance in function value, quite useful'
    , tolfunhist='1e-12  #v termination criterion: tolerance in function value history'
    , tolstagnation='int(100 + 100 * N**1.5 / popsize)  #v termination if no improvement over tolstagnation iterations'
    , tolx='1e-11  #v termination criterion: tolerance in x-changes'
    , transformation='None  # [t0, t1] are two mappings, t0 transforms solutions from CMA-representation to f-representation (tf_pheno), t1 is the (optional) back transformation, see class GenoPheno'
    , typical_x='None  # used with scaling_of_variables'
    , updatecovwait='None  #v number of iterations without distribution update, name is subject to future changes' # TODO: rename: iterwaitupdatedistribution?
    , verb_append='0  # initial evaluation counter, if append, do not overwrite output files'
    , verb_disp='100  #v verbosity: display console output every verb_disp iteration'
    , verb_filenameprefix='outcmaes  # output filenames prefix'
    , verb_log='1  #v verbosity: write data to files every verb_log iteration, writing can be time critical on fast to evaluate functions'
    , verb_plot='1  #v in fmin(): plot() is called every verb_plot iteration'
    , verb_time='True  #v output timings on console'
    , vv='0  #? versatile variable for hacking purposes, value found in self.opts[\'vv\']'
     ):

    try: # pass on KeyboardInterrupt
        opts = locals()  # collect all local variables (i.e. arguments) in a dictionary
        del opts['func'] # remove those without a default value
        del opts['args']
        del opts['x0']      # is not optional, no default available
        del opts['sigma0']  # is not optional for the constructor CMAEvolutionStrategy
        if not func:  # return available options in a dictionary
            return Options(opts, True)  # these opts are by definition valid

        # TODO: this is very ugly:
        incpopsize = Options({'incpopsize':incpopsize}).eval('incpopsize')
        restarts = Options({'restarts':restarts}).eval('restarts')
        del opts['restarts']
        noise_handling = Options({'noise_handling': noise_handling}).eval('noise_handling')
        del opts['noise_handling']# otherwise CMA throws an error

        irun = 0
        best = BestSolution()
        while 1:
            # recover from a CMA object
            if irun == 0 and isinstance(x0, CMAEvolutionStrategy):
                es = x0
                x0 = es.inputargs['x0']  # for the next restarts
                if sigma0 is None or not np.isscalar(array(sigma0)):
                    sigma0 = es.inputargs['sigma0']  # for the next restarts
                # ignore further input args and keep original options
            else:  # default case
                if irun and opts['restart_from_best']:
                    print('CAVE: restart_from_best is typically not useful')
                    es = CMAEvolutionStrategy(best.x, sigma0, opts)
                else:
                    es = CMAEvolutionStrategy(x0, sigma0, opts)
                if opts['eval_initial_x']:
                    x = es.gp.pheno(es.mean, bounds=es.gp.bounds)
                    es.best.update([x], None, [func(x, *args)], 1)
                    es.countevals += 1

            opts = es.opts  # processed options, unambiguous

            append = opts['verb_append'] or es.countiter > 0 or irun > 0
            logger = CMADataLogger(opts['verb_filenameprefix'], opts['verb_log'])
            logger.register(es, append).add()  # initial values, not fitness values


            noisehandler = NoiseHandler(es.N, noise_handling, np.median, opts['noise_reevals'], opts['noise_eps'], opts['eval_parallel'])
            while not es.stop():
                X, fit = es.ask_and_eval(func, args, evaluations=noisehandler.evaluations,
                                         aggregation=np.median) # treats NaN with resampling
                # TODO: check args and in case use args=(noisehandler.evaluations, )

                if 11 < 3 and opts['vv']:  # inject a solution
                    # use option check_point = [0]
                    if 0 * np.random.randn() >= 0:
                        X[0] = 0 + opts['vv'] * es.sigma**0 * np.random.randn(es.N)
                        fit[0] = func(X[0], *args)
                        # print fit[0]
                es.tell(X, fit)  # prepare for next iteration
                if noise_handling:
                    es.sigma *= noisehandler(X, fit, func, es.ask, args)**opts['noise_change_sigma']
                    es.countevals += noisehandler.evaluations_just_done  # TODO: this is a hack, not important though

                es.disp()
                logger.add(more_data=[noisehandler.evaluations, 10**noisehandler.noiseS] if noise_handling else [],
                           modulo=1 if es.stop() and logger.modulo else None)
                if opts['verb_log'] and opts['verb_plot'] and \
                    (es.countiter % max(opts['verb_plot'], opts['verb_log']) == 0 or es.stop()):
                    logger.plot(324, fontsize=10)

            # end while not es.stop
            mean_pheno = es.gp.pheno(es.mean, bounds=es.gp.bounds)
            fmean = func(mean_pheno, *args)
            es.countevals += 1

            es.best.update([mean_pheno], None, [fmean], es.countevals)
            best.update(es.best)  # in restarted case

            # final message
            if opts['verb_disp']:
                srestarts = (' after %i restart' + ('s' if irun > 1 else '')) % irun if irun else ''
                for k, v in list(es.stop().items()):
                    print('termination on %s=%s%s (%s)' % (k, str(v), srestarts, time.asctime()))

                print('final/bestever f-value = %e %e' % (es.best.last.f, best.f))
                if es.N < 9:
                    print('mean solution: ' + str(es.gp.pheno(es.mean)))
                    print('std deviation: ' + str(es.sigma * sqrt(es.dC) * es.gp.scales))
                else:
                    print('mean solution: %s ...]' % (str(es.gp.pheno(es.mean)[:8])[:-1]))
                    print('std deviations: %s ...]' % (str((es.sigma * sqrt(es.dC) * es.gp.scales)[:8])[:-1]))

            irun += 1
            if irun > restarts or 'ftarget' in es.stopdict or 'maxfunevals' in es.stopdict:
                break
            opts['verb_append'] = es.countevals
            opts['popsize'] = incpopsize * es.sp.popsize # TODO: use rather options?
            opts['seed'] += 1

        # while irun

        es.out['best'] = best  # TODO: this is a rather suboptimal type for inspection in the shell
        if 1 < 3:
            return es.result() + (es.stop(), es, logger)

        else: # previously: to be removed
            return (best.x.copy(), best.f, es.countevals,
                    dict((('stopdict', CMAStopDict(es.stopdict))
                          ,('mean', es.gp.pheno(es.mean))
                          ,('std', es.sigma * sqrt(es.dC) * es.gp.scales)
                          ,('out', es.out)
                          ,('opts', es.opts)  # last state of options
                          ,('cma', es)
                          ,('inputargs', es.inputargs)
                          ))
                   )
        # TODO refine output, can #args be flexible?
        # is this well usable as it is now?
    except KeyboardInterrupt:  # Exception, e:
        if opts['verb_disp'] > 0:
            print(' in/outcomment ``raise`` in last line of cma.fmin to prevent/restore KeyboardInterrupt exception')
        raise  # cave: swallowing this exception can silently mess up experiments, if ctrl-C is hit

Mh = Misc.MathHelperFunctions

def run_bbob_bench(fun, dim, ftarget, maxfunevals):
#    fmin(fun, np.ones(dim)*0.1, 0.1)
    optim = CMAEvolutionStrategy(np.ones(dim)*0.1, 0.1)
#    optim.optimize(fun, ftarget, maxfunevals)
    optim.optimize_ViE(fun, ftarget, maxfunevals)


class CMAVIE:
    
    def __init__(self):
        pass
        
    def add_keywords(self,params):
        
        # CMAViE calculation flags
        params.add('maxFevals','maxFevals','int',16000)
        params.add('evalByRestarts','evalByRestarts','int',4000)
        params.add('continuousRestarts', 'continuousRestarts', 'str', 'off')
        params.add('updateSimultaneous', 'updateSimultaneous', 'str', 'off')
        params.add('writeConverged', 'writeConverged', 'str', 'off')
        
        return params
    
    #just a function to check the optimizer-specific keywords passed via the input file
    def check_keywords(self):
        
        # For CMAVIE only, decide whether we want to restart from previous best
        if self.params.continuousRestarts == 'on' :
            self.params.continuousRestarts=True
        elif self.params.continuousRestarts== 'off':
            self.params.continuousRestarts=False
        else:
            print ('ERROR: continuousRestarts should either be "on" or "off"!')
            sys.exit(1)   
            
        # check the constraints which makes you write file in CMAViE only if boundary constraints are satisfied
        if self.params.writeConverged == 'on' :
            self.params.writeConverged=True
        elif self.params.writeConverged== 'off':
            self.params.writeConverged=False
        else:
            print ('ERROR: constraint should either be "on" or "off"!')
            sys.exit(1) 

        # Simultaneous constraint optimization flag
        if self.params.updateSimultaneous == 'on':
            self.params.updateSimultaneous = True
        else:
            self.params.updateSimultaneous = False
        if self.params.updateSimultaneous and not self.constrained_opt:
            print ('ERROR: simultaneous boundary optimization was switched, though the constraint optimization flag was off...')
            sys.exit(1)


    def launch(self,params,space,fitness):
        
        # Equivalent of __init__(self)
        self.params = params
        self.space = space
        self.fitness = fitness

        self.func = fitness.evaluate
        
        self.space.cell_size=self.space.high-self.space.low
        self.params.dimensions=len(self.space.cell_size)
        # ------------------------- #
        
        '''standard function used by POW to launch the optimizer'''
        if rank == 0:
            import subprocess
            subprocess.call('rm CMAlog.txt',shell=True) 
            subprocess.call('rm log.txt',shell=True) 
            
            #master node checks for optimizer-specific keywords consistency
            self.check_keywords()

        # for now this part is off because we want to use multiple restarts for the code
        if not self.params.constrained_opt and 11 < 3:
            self.optimizer = CMAEvolutionStrategy(np.zeros(len(self.space.low))+0.5, 0.1, {"tolx":1e-4})
            self.optimizer.optimize_ViE_POW(self.func, None, self.params, self.space)
            
        else:            
            function_evals_counter = 0
            max_function_evals = self.params.maxFevals
            bool = 0

            if rank == 0:

                # keep track of the success
                success_counter = 0
                restarts = 0
                # __NEW__
#                 self.params.continuousRestarts=True
#                     
#             else:
#                 self.params.continuousRestarts=None
                
#             comm.Barrier()
#             self.params.continuousRestarts = comm.bcast(self.params.continuousRestarts, root = 0)
#             comm.Barrier()

            while function_evals_counter < max_function_evals:
                
                if not self.params.continuousRestarts:
                    # creating instance of optimizer with mean randomly distibuted around 0 and 1.
                    self.optimizer = CMAEvolutionStrategy(np.random.uniform( 0.0, 1.0, len(self.space.low)), 0.1, {"tolx":1e-3}) 
                
                else:
                    # at the first run there is no loaded restart coordinates
                    if function_evals_counter == 0:
                        self.optimizer = CMAEvolutionStrategy(np.random.uniform( 0.0, 1.0, len(self.space.low)), 0.1, {"tolx":1e-3}) 
                    elif bool == 0:
                        self.optimizer = CMAEvolutionStrategy(np.random.uniform( 0.0, 1.0, len(self.space.low)), 0.1, {"tolx":1e-3})
                    else:
                        if rank == 0:
                            # open the params file and read the last input
                            mat=np.loadtxt("CMAVIE_coords.dat")
                            last=mat[0]
                            
                        else:
                            last=None
                        
                        comm.Barrier()
                        last = comm.bcast(last, root = 0)
                        comm.Barrier()
                        
                        self.optimizer = CMAEvolutionStrategy(last, 0.05, {"tolx":1e-3}) 
                
                
                if 1 < 3:
                    # POW with CMAVIE optimization
                    bool, feval, X, fitvals = self.optimizer.optimize_ViE_POW(self.func, None, self.params, self.space)
                else:
                    # POW with simple CMAES optimization
                    self.optimizer.optimize_CMA_POW(self.func, None, self.params, self.space)
                    feval = self.params.maxFevals
                    bool=1

                if rank == 0:
                    # increase the fevals counter:
                    function_evals_counter += feval
                    restarts += 1

                    if bool == 1:
                        success_counter += 1

                comm.Barrier()
                function_evals_counter = comm.bcast(function_evals_counter, root = 0)
                comm.Barrier()

            if rank == 0:
                print (">>>> TOTAL NUMBER OF FUNCTION EVALS: "+str(function_evals_counter))
                print (">>>> successes/restarts: "+str(success_counter)+" / "+str(restarts))
                
#                 if self.params.continuousRestarts:
#                     self.params.CMAVIE_coords.close()


if __name__ == "__main__":
#
    f = FitnessFunctions()
    func = f.rastrigin

    # import the bbob functions:
#    function_ids = bbobbenchmarks.nfreeIDs
#    instance = 1
#    opts = dict(algid='', comments='no comments')
#    datapath = ''
#    functionNo = 0
#
#    iinstance = 1
#    f = fgeneric.LoggingFunction(datapath, **opts)
#    f.setfun(*bbobbenchmarks.instantiate(function_ids[functionNo], iinstance=iinstance))
##
###    fmin(f.rastrigin, np.ones(20)*0.1, 0.1)
###    func = f.hyperelli
#    func = f.evalfun


    result1 = []
    result2 = []
    for i in xrange (0,1,1):
        print(i)

        optim = CMAEvolutionStrategy(np.ones(5)*0.1, 0.1)
        result1.append(optim.optimize_ViE(func, None, None))
#        print("\n")
#        print ("-----------------------------------------------------------------------------------------------")
#        print("\n")
        # other way to optimise:
        optim1 = CMAEvolutionStrategy(np.ones(5)*0.1, 0.1)
        result2.append(optim1.optimize(func))

#        print ("\n\n")

    print ("\n\n===================================================================================================\n\n")
#    print ("array cma vie: "+str(result1))
#    print ("\narray cma norm: "+str(result2))
    print ("\n>> cma vie:\t"+str(sum(result1)))
    print (">> cma norm:\t"+str(sum(result2)))

   # fmin(func, np.ones(20)*0.1, 0.1)
