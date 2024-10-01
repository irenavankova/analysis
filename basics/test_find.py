#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate random 3D data points
x = np.linspace(11,20,10)
y = x/0
y0 = x*0
y1 = x*1
print(x)
ind_x = np.where(x > 15)
print(ind_x)
bool_x = x > 15
print(bool_x)

print('INDEX WHERE')
ind_x = np.where((x > 15) & (x < 19))
print(ind_x)
print(x[ind_x])

print('BOOL WHERE')
bool_x = (x > 15) & (x < 19)
print(bool_x)
print(x[bool_x])

print('INDEX CONDITION')
whr_x = np.where((x > 15) & (x < 19), x, y)
print(whr_x)

whr_x = np.where((x > 15) & (x < 19), y, x)
print(whr_x)

ind_x = np.where((x > 15) & (x < 19), 1, 0)
print(ind_x)
print(x*ind_x)




