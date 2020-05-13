import numpy as np

def constraint_check(multimer):


    width = multimer.get_width()
    if (width > 58 and width < 62):
       width=60

    height = multimer.get_height()
    if (height > 45 and height < 49):
       height=47

    #select atoms
    ca1 = multimer.atomselect(2, "A", 10, "OD2")
    ca2 = multimer.atomselect(3,"A",105,"NE2")
    ca3 = multimer.atomselect(2, "A", 64, "NZ")
    ca4 = multimer.atomselect(3,"A",55,"NH1")

    dist1 = multimer.distance(ca1, ca2)
    if dist1>5 and dist1<13:
        dist1=9
    
    dist2 = multimer.distance(ca3, ca4)
    if dist2>6 and dist2<14:
        dist2=10

    #return dist1, dist2
    return width, height, dist1, dist2
