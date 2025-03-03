#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate random 3D data points

A = np.random.randint(0, 100, size=(5, 3))

print(A.shape)
print(A)

#bool_x = (A > 0.5)
#print(bool_x)

a = [0,1]
b = [1,2]
B = A[a,b]

print(B)

iis = [True, False, False, True]
#iis = (A > 50)
print(iis)
ix = np.where(iis)[0]
print(ix)








