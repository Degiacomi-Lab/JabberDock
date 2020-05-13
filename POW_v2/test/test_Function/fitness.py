import numpy as np

class Fitness:
    def __init__(self,data,params):
        pass

    def evaluate(self,num,pos):
        #Rastrigin function
        return 10*len(pos)+np.sum(pos**2-10*np.cos(2*np.pi*pos))


