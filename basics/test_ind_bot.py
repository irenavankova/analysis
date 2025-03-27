#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate random 3D data points

A = np.array([[[1, 2], [3, 4]],[[5, 6], [7, 8]]])

print('A')
print(A.shape)
print(A)

B = np.array([[0, 1],[0,1]])

#rows, cols = np.indices(B.shape)  # Create row and column indices for B
#result = A[rows, cols, B]
result = A[np.arange(B.shape[0])[:, None], np.arange(B.shape[1]), B]
print('A(B)')
print(result.shape)
print(result)








