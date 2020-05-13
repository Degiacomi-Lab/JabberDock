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


from Default import Parser as P
from Default import Space as S
from Default import Postprocess as PP

import numpy as np
import sys

class Parser(P):
    def __init__(self):
        pass
    def check_variables(self):
        pass

class Data:
    def __init__(self,params):
        pass      

class Space(S):
    def __init__(self,params,data):

        #assign low boundaries
        if params.low_input!="NA" :
            self.low=np.zeros(len(params.low_input))
            for i in xrange(0,len(params.low_input),1):
                self.low[i]=params.low_input[i]
        else:
            print "ERROR: boundaryMin should be defined"
            sys.exit(1) 
        
        #assign high boundaries
        if params.high_input!="NA" :
            self.high=np.zeros(len(params.low_input))
            for i in xrange(0,len(params.high_input),1):
                self.high[i]=params.high_input[i]
        else:
            print "ERROR: boundaryMax should be defined"
            sys.exit(1)
                              
        #set boundary type (default is periodic)
        self.boundary_type=np.zeros(len(params.low_input))
        if params.boundary_type!="NA":
            for i in xrange(0,len(params.low_input),1):
             	self.boundary_type[i]=params.boundary_type[i]

class Fitness:
    def __init__(self,data,params):
        pass
    def evaluate(self,num,pos):
        #rastrigin function (used for benchmark)
        return 10*len(pos)+np.sum(pos**2-10*np.cos(2*np.pi*pos))


class Postprocess(PP):

    def __init__(self,data,params):
        pass

    def run(self) :
        pass
