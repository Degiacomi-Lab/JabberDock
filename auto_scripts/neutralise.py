#!/usr/bin/env python

from __future__ import division

""" A short script for counting and neutralising a simulation setup"""

import numpy as np
import sys 

filename = sys.argv[1]

with open(filename, "r") as ins:
    array = []
    for line in ins:
        array.append(float(line))

charge = np.sum(array)

print(int(charge))
