# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 10:48:46 2014

@authors:   Giorgio Tamo - giorgiotamo@gmail.com
            Alexanderbelushkin - alexander.belushkin@epfl.ch
"""


#------------IMPORTANT!---------------------------
#---All attribute arrays in the input-------------
#---objects, i.e. vFunction, problemOpt, inopts---
#---have to be of the numpy.array type------------
#------------IMPORTANT!---------------------------

from copy import deepcopy as decc
import numpy
import random
from numpy import array as nua
from numpy import inf
import math

def algo_de_epsilon(vFunction, problemOpt, inopts):
    # Init DE parameters
    class params_init:
        def __init__(self, vFunction, problemOpt, inopts):
            self.popSize = inopts.PopSize
            self.maxGen = inopts.MaxGen
            self.LocalCount = inopts.LocalCount#Max local fevals
            self.CR = inopts.CR
            self.F = inopts.F
            self.F0 = self.F
            self.stopFitness = inopts.StopFitness
            self.select_de_best = False
            self.dim = problemOpt.n
            self.binaryCrossover = inopts.binaryCrossover
            self.lb = problemOpt.domain[:, 0]#1
            self.ub = problemOpt.domain[:, 1]#2
            self.epsilon0 = 0
            self.epsilonTheta = 0.2
            self.epsilonLambda = 10^-5
            self.epsilonTc = 1000000 # _GT_ was 2500
            self.epsilonTlambda = 0.95*self.epsilonTc
            self.epsilon = 0
            self.ne = 3
            self.pg = 0.01
            self.rg = 3
            self.cp = 5
            self.vFunction = vFunction
            self.tolx = 10e-6
            self.constraint_target = inopts.target
    
    p_local=params_init(vFunction, problemOpt, inopts)
    global params #make params global so that don't have to pass to each fuction
    params=p_local
    
    # Init constraint parameters
    target = numpy.array( [[-inf, -inf]] )
    
    target=numpy.concatenate((target, problemOpt.target[1:]),axis=0)
    
    
    bcIndexesLb = numpy.where(problemOpt.domain[:, 0] != -inf)#1
    bcIndexesUb = numpy.where(problemOpt.domain[:, 1] != inf)#2
    
    global constraintEvalCount #should be an object!----------------------line 41 
    global objEvalCount#should be an object!
    constraintEvalCount = 0
    objEvalCount = 0

    # Init best individual
    class best:#--------------------------------------line 47
        def __init__(self):
            self.fitness = inf
            self.x = []
            self.ix = 0
            self.error = None
            self.fEvals=None
            self.cEvals=None
            
    class pop_member():
        def __init__(self):
            self.x=None
            self.fitness=None
            self.ix=None
            self.fEvals=None
            self.cEvals=None
            
    class trial():
        def __init__(self):
            self.x=None
            self.fitness=None
            self.ix=None
            self.fEvals=None
            self.cEvals=None

    class evalCount():
        def __init__(self):
            self.fevals=None
            self.cavals=None
        
    evalCount = evalCount()
    best=best()
    pop_m=pop_member()        
     
    population=[]#---------------line52
    for i in range(params.popSize):
        population.append(pop_member())    
    
       
    # Generate initial population and evaluate it
  
    for i in range(params.popSize):
          
        # _GT_ generate individual
        difference = [ub-lb for ub,lb in zip(inopts.ub, inopts.lb)]
        randm = numpy.random.random( problemOpt.n ) #[problemOpt.n,1]         
        transient_array = nua([a*b for a,b in zip(difference,randm) ])
        population[i].x = nua([a+b for a,b in zip(inopts.lb,transient_array) ])          
        
        # _GT_ evaluate individual
        constraints, fitness = vFunction(0, population[i].x)
        
        # -----POW  DE implementation in constrained optimization
        # -> POW constraints and objective function have to be transformed
        
        result = [];
        result.append(constraints[-1])
    
        # transform parameters into constraint violations
        for t in xrange(len(params.constraint_target)):
            c  = constraints[t]
            lb = params.constraint_target[t]
            ub = params.constraint_target[t]

            lerr = lb - c
            result.append(lerr)
            uerr = c - ub
            result.append(uerr)  
        # ----------
                
        vValues = nua(result)
        vValues = vValues.T        
        
        # _GT_ compute constraints violations
        cviol = nua( vValues[1:] )
        
        problemOpt.domain=nua(problemOpt.domain)
        bcIndexesLb=nua(bcIndexesLb)
        bcIndexesUb=nua(bcIndexesUb)
        bviol_first_row = nua( [ problemOpt.domain[bcIndexesLb,0] - population[i].x[bcIndexesLb] ] )
        bviol_second_row = nua( [ population[i].x[bcIndexesUb] - problemOpt.domain[bcIndexesUb, 1] ] )
        bviol = nua( [ bviol_first_row, bviol_second_row ] )
          
        # _GT_ update params of individual
          
        population[i].fitness = vValues[0]
        population[i].error = sum( cviol[cviol>0] ) + sum( bviol[bviol>0] )         
        population[i].ix = i        
        population[i].fEvals = False
        population[i].cEvals = False
        
        #--------------------------------------------------------------------
        # _GT_ extract best of initial population
        #    if result of epsilon comparison is 1 individual is better than
        #    best and becomes new best

        if i == 0:
           res = 1
        else:
                      
           [res, population[i], best ] = Better(population[i], best)
        

        if res==1:
           best.x = population[i].x
           best.fitness = population[i].fitness
           best.ix = i
           best.error = population[i].error
           best.fEvals = population[i].fEvals
           best.cEvals = population[i].cEvals
  
    # Epsilon control - Eps0 - _GT_ extract all errors and update Epsilon0
    
    errors = numpy.zeros(params.popSize)
    
    
    for i in range(params.popSize):
       errors[i] = population[i].error#----------------------------------------line99
    
    #[trash, errorsIx] = sort(errors) # _GT_ extract sorted indices of errors

    errorsIx=sorted(range(len(errors)),key=errors.__getitem__)

    # params.epsilon0 = errors(errorsIx(round(params.epsilonTheta*params.popSize)))
    params.epsilon0 = errors[errorsIx[0]] * params.epsilonTheta   

    # Elite preservation
    elites = [[]]*params.ne # _GT_ ne is number of elite individuals
    elites[0] = decc(population[errorsIx[0]])
    elites[1] = decc(population[errorsIx[1]])
    elites[2] = decc(population[errorsIx[2]])

    
    # _GT_ =========== MAIN LOOP ===========
#     Dist_conv=[[],[]]    
    for g in range(int(params.maxGen)):
         
        # _GT_ stopping criteria ------------------------------
           
        
        if (best.error == 0) and (best.fitness < params.stopFitness):
            xmin = best.x
            vMin = best.fitness
            evalCount.fevals = objEvalCount
            evalCount.cevals = constraintEvalCount
            stopflag = [ 'fitness' ]
            return [xmin, vMin, evalCount, stopflag, Dist_conv]
        
        # -------------------------------------------------
    
        if objEvalCount>params.LocalCount:
            xmin = best.x
            vMin = best.fitness
            evalCount.fevals = objEvalCount
            evalCount.cevals = constraintEvalCount
            stopflag = [ 'Local counts' ]
            return [xmin, vMin, evalCount, stopflag, Dist_conv]
            
        # Epsilon control - EpsT
        if g < params.epsilonTc:
            
#             if g > params.epsilonTlambda:
#                 params.F = 0.3*params.F0 + 0.7
#                 params.cp = 0.3*params.cp + 0.7*3
#             else:
#                 params.cp = (-5 - math.log(params.epsilon0) )/math.log(0.05)
#             if params.cp < 3:
#                params.cp = 3
#             
#             params.cp = 5
#             params.epsilon = params.epsilon0*(1- (g/params.epsilonTc))**params.cp
            params.epsilon = 0
        else:
            params.epsilon = 0
            # Disable elites
            elites = []
            params.ne = 0

        for i in range(params.popSize):
            # store old x value to make sure not to stall
#             if trial != None:
#                 oldbestX = trial.x
#                 oldbestF = trial.fitness
            
                 
            [p1, p2, p3] = de_select(i, params, best)

            if params.ne > 0:
                popForRecomb = numpy.concatenate( (population, elites) )
            else:
                popForRecomb = decc(population)

            if params.binaryCrossover:
                   
                trial.x = de_1_Bin(population[i].x, popForRecomb[p1].x, popForRecomb[p2].x, popForRecomb[p3].x, params)
            else:
                trial.x = de_1_Exp(population[i].x, popForRecomb[p1].x, popForRecomb[p2].x, popForRecomb[p3].x, params)
                
            # Keep solution in box constraints
            trial.x = rangeKeeperReflection(trial.x, params)
                        
            trial.x = trial.x.T
            
            # Evaluate TRIAL individual
            #--------------------------------------------------------------------
            constraints, fitness = vFunction(0, trial.x)
        
            # -----POW  DE implementation in constrained optimization
            # -> POW constraints and objective function have to be transformed
            
            result = [];
            result.append(constraints[-1])
        
            # transform parameters into constraint violations
            for t in xrange(len(params.constraint_target)):
                c  = constraints[t]
                lb = params.constraint_target[t]
                ub = params.constraint_target[t]
    
                lerr = lb - c
                result.append(lerr)
                uerr = c - ub
                result.append(uerr)  
            # ----------
                    
            vValues = nua(result) 
              
            vValues = vValues.T
            
            cviol = vValues[1:]
            
            bviol_first_row = nua( [ problemOpt.domain[bcIndexesLb,0] - trial.x[bcIndexesLb] ] )
            
            bviol_second_row = nua( [ trial.x[bcIndexesUb] - problemOpt.domain[bcIndexesUb, 1] ] )
            
            
            bviol = nua( [ bviol_first_row, bviol_second_row ] )
            
          
            trial.error = sum(cviol[cviol > 0]) + sum(bviol[bviol > 0])
                      
            trial.fitness = vValues[0]
            trial.fEvals = False
            trial.cEvals = False
          
                      
            [res, trial, population[i]] = Better(trial, population[i])
                        
            
            if res==1:
                population[i].fitness = trial.fitness
                population[i].x  = trial.x
                population[i].error = trial.error
                population[i].fEvals = trial.fEvals
                population[i].cEvals = trial.cEvals
          
          
            # Update best solution
                   

        
            [res, trial, best] = Better(trial, best)
                    
                        
            
            if res==1:
                best.x = trial.x
                best.ix = i
                best.fitness = trial.fitness
                best.error = trial.error
                best.fEvals = trial.fEvals
                best.cEvals = trial.cEvals
                
                
            # Update elites
            elites.sort(BetterForSort)#params are global
           
            # Check if trial best than worst elite
            [res, trial, elites[params.ne-1]] = Better(trial, elites[params.ne-1])

            if res==1:
                elites[params.ne-1].fitness = trial.fitness
                elites[params.ne-1].x  = trial.x
                elites[params.ne-1].error = trial.error
                elites[params.ne-1].fEvals = trial.fEvals
                elites[params.ne-1].cEvals = trial.cEvals

            # --------- outprint
            print "> Best - Fitness : "+str(best.fitness)+" - Error : "+str(best.error)
            # evaluating whether stalling termination criteria has been reached
#             dist = numpy.linalg.norm(oldbestX-best.x)
#             print dist
#                 print (oldbestX - best.x)
#                 if numpy.fabs(oldbestX - best.x) < params.tolx:
#                     print "TOLX!!!!"
            
                
            
            

    # _GT_ stopping criteria ------------------------------

            
    xmin = best.x
    vMin = best.fitness
    evalCount.fevals = objEvalCount
    evalCount.cevals = constraintEvalCount
    stopflag = [ 'MaxEvals' ]   
    return [xmin, vMin, evalCount, stopflag]
    #END OF algo_de_epsilon function
          
#      =================================================================
# _GT_ =========================== FUNCTIONS ===========================
#      =================================================================

def setdiff(a,b):
    res=[]    
    for i in a:
        if i not in b:
            res.append(i)
    res.sort()
    return res        

def de_select(i, params, best):
     if params.ne > 0:
         possibleIx = range(params.popSize+params.ne)
     else:
         possibleIx = range(params.popSize)
         
     possibleIx = setdiff(possibleIx, [i])
        

     if params.select_de_best:
         p1 = decc(best.ix)
     else:
         p1 = possibleIx[random.randint(0,len(possibleIx)-1)]
    
     possibleIx = setdiff(possibleIx, [p1])
     p2 = possibleIx[random.randint(0,len(possibleIx)-1)]
    
     possibleIx = setdiff(possibleIx, [p2])
     p3 = possibleIx[random.randint(0,len(possibleIx)-1)]
     return [p1, p2, p3]

def de_1_Bin(xi, x1, x2, x3, params):
    # original
    #             j = randi(params.dim, 1)
    #             for l=0:dim-1
    #                if ( l == j || rand(1) < params.CR)
    #                   trial(l) = x1(l) + params.F*(x2(l)-x3(l))
    #                else
    #                   trial(l) = xi(l)
    #                end
    #             end

    # vectorized
    j = numpy.random.randint(0,params.dim-1)
    r = numpy.random.rand(params.dim-1)
    
    ix1 = r < params.CR#boolean array, True and False values
    ix2 = ~ix1
    
    # the order is important here!
    trial=numpy.zeros(len(ix1))   
    
    for i in range(len(ix2)): 
        if ix2[i]:
            trial[i]=decc(xi[i])
            
    for i in range(len(ix1)):
        if ix1[i]:
            trial[i] = x1[i] + params.F*(x2[i]-x3[i])
    
    trial[j] = x1[j] + params.F*(x2[j]-x3[j])
     
    return trial
         
#x_child generation
def de_1_Exp(xi,x1,x2,x3,params):
   
    j=numpy.random.randint(0,params.dim)
    l=0
    # First instance of do-while
    trial=numpy.zeros(len(x1))
    
    
    trial[j] = x1[j] + params.F*(x2[j]-x3[j])       
    j = j + 1
    if j>params.dim-1:
        j=0
    l = l + 1
    while (l<params.dim) and (numpy.random.rand() < params.CR):
        trial[j] = x1[j] + params.F*(x2[j]-x3[j])
        j = j + 1
        if j > params.dim-1:
            j = 0           
        l = l + 1
       
    while l < params.dim:
        trial[j] = decc(xi[j])
        j = j + 1
        if j > params.dim-1:
            j = 0
        l = l + 1
    
    
    return trial        
   
def rangeKeeperReflection(xx, params):
    x=xx    
    
    for j in range(len(x)):
        if x[j]<params.lb[j]:
            exceed = decc(params.lb[j] - x[j])
            if exceed > (params.ub[j]-params.lb[j]):
                exceed = exceed - round(exceed/(params.ub[j]-params.lb[j]))*(params.ub[j]-params.lb[j])
            x[j] = decc(params.lb[j]+exceed)
        elif x[j]>params.ub[j]:
            exceed =  decc(x[j] - params.ub[j])
            if exceed > (params.ub[j]-params.lb[j]):
                exceed = exceed - round(exceed/(params.ub[j]-params.lb[j]))*(params.ub[j]-params.lb[j])
            x[j] = decc(params.ub[j] - exceed)
         
    return x


def Better(p1,p2):
    res, p11, p22 = epsilonComparison(p1,p2)
    return [res, p11, p22]


def BetterForSort(p1, p2):
    [res, trash1, trash2] = Better(p1, p2)
    res=res*-1    
    return res

def epsilonComparison(p1, p2):
    global constraintEvalCount
    global objEvalCount
   
    global params


    p11=decc(p1)
    p22=decc(p2)
    
    
    
    epsilon=params.epsilon
    
    e1 = decc(p11.error)
    e2 = decc(p22.error)
    
    if not p11.cEvals:
        p11.cEvals = True # _GT_ switching c evaluation bool      
        constraintEvalCount = constraintEvalCount + 1
    
    if not p22.cEvals:
        p22.cEvals = True # _GT_ switching c evaluation bool      
        constraintEvalCount = constraintEvalCount + 1
        
    if (epsilon == 0) or (e1 > epsilon) or (e2 > epsilon):
       if e1 < e2:
          res = True
       elif e1 > e2:
          res = False
       else:
          # _GT_ update individual
          if  not p11.fEvals:
             p11.fEvals = True
             objEvalCount = objEvalCount + 1
          
          if not p22.fEvals:
             p22.fEvals = True
             objEvalCount = objEvalCount + 1
          
          
          if p11.fitness < p22.fitness:
             res = True
          else:
             res = False
          
    else:
       
       # _GT_ update individuals
       if not p11.fEvals:
          p11.fEvals = True
          objEvalCount = objEvalCount + 1
       
       if not p22.fEvals:
          p22.fEvals = True
          objEvalCount = objEvalCount + 1
       

       if p11.fitness != p22.fitness:
          if p11.fitness < p22.fitness:
             res = True
          else:
             res = False
          
       else:
          if e1 < e2:
              res = True
          else:
              res=False
    if res:
        res = 1
    else:
        res = -1
   
    
    return [res, p11, p22]              
      
#END_OF_PROGRAM

# ======================== LAUNCHING FROM POW ========================

class DE_epsilon:
    def __init__(self):
        pass
    
    def add_keywords(self,params):
        params.add('maxFevals','maxFevals','int',16000)
        params.add('repeat','repeat','int',1)
        
        return params
    
    def check_keywords(self):
        pass
    
    def launch(self,params,space,fitness):
        ''' Standard POW function to launch optimizer '''
        
        
        # Equivalent of __init__(self)
        self.params = params
        self.space = space
        self.fitness = fitness

        self.func = fitness.evaluate
        
        self.space.cell_size=self.space.high-self.space.low
        
         
        self.domain = numpy.array([[self.space.low[i],self.space.high[i]] for i in xrange(0,len(self.space.high),1) ])
        
        self.params.dimensions=len(self.space.cell_size)
        
        self.target = [[-inf, 0]]*len(self.params.target)
        self.target.insert(0,-inf)
        
        class fParams():
            def __init__(self, id, Nequalities, Ninequalities):
                
                self.fId = id
                fParams.gn = Nequalities 
                fParams.hn = Ninequalities
        self.fParams=fParams('DE_for_POW', len(self.target), 0)
        # ------------------------- #
        
        
        # Set the problem options
        class problemOpt():
            def __init__(self,n, domain, initDomain, target, fopt, startPoint, fParams):
                
                self.n          = n
                self.domain     = domain 
                self.initDomain = initDomain 
                self.target     = target 
                self.fopt       = fopt 
                self.startPoint = startPoint 
                self.fParams    = fParams 
#                 self.xopt       = xopt 
                
        
        
        # give POW parameters to DE
        problemOpt = problemOpt(self.params.dimensions, self.domain, self.domain, self.target, -numpy.inf, [], self.fParams) 
        
          
        # Set algorithm options
        class inopts():
            def __init__(self, domain, maxFevals, target):
                        
                self.StopFitness    = -numpy.inf
                self.MaxFunEvals    = maxFevals
                self.lb             = domain[:, 0] 
                self.ub             = domain[:, 1] 
                self.PopSize        = 15 
                self.CR             = 0.9 
                self.F              = 0.7 
                self.MaxGen         = numpy.floor(self.MaxFunEvals/self.PopSize) 
                self.LocalCount     = 20000 #Maximum number of fevals in one call of DE
                self.binaryCrossover= False 
                self.target         = target
        
        inopts = inopts(self.domain, self.params.maxFevals, self.params.target )      
        # Reset objective and constraints function evaluations
        nTrials = self.params.repeat
        fevals = numpy.zeros( nTrials) 
        cevals = numpy.zeros(nTrials) 
        startevals = numpy.zeros( nTrials) 
        success = numpy.zeros(nTrials) 
        
        #Doing just 1 trial
        
        pOpt = decc(problemOpt) 
        iOpts = decc(inopts) 
              
        # Until fEvals left and problem not solved, re-iterate
        fcount = 0 
        ccount = 0 
        curRestart = 0 
        
#         from DE_epsilon import algo_de_epsilon

        #[ xMin, vMin, evalCount, stopFlag ] = algo_de_epsilon( bFunction, pOpt, iOpts) 
        
        # wrt=open(curBench,'w') 
        while iOpts.MaxGen > 0:
            [ xMin, vMin, evalCount, stopFlag] = algo_de_epsilon( self.func, pOpt, iOpts) 
            fcount += evalCount.fevals 
            ccount += evalCount.cevals 
            iOpts.MaxGen -= numpy.ceil(evalCount.cevals/iOpts.PopSize)   
            if stopFlag[0] == 'fitness':
                break
            curRestart = curRestart + 1 
        #     wrt.write(str(Distance_conv[0])) 
        #     wrt.write('\n') 
        #     wrt.write(str(Distance_conv[1])) 
            
        #     wrt.close() 
                    
        print 'xMin     = '+str(xMin) 
        print 'vMin     = '+str(vMin) 
        print 'cEvals   = '+str(evalCount.cevals)
        print 'fEvals   = '+str(evalCount.fevals)
        print 'stopFlag = '+str(stopFlag)