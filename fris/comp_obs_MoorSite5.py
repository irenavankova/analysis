#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

# 1. Load the .mat file
mat_data = sio.loadmat("/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/DOVuFRIS/Fris_Apres/code/ocean_data/Site5_ocean_data/Mooring5c/s5cInterp.mat")
s5c = mat_data["s5c"].squeeze()

# 2. Define indices (MATLAB: 1,3,4,5,7,8) -> Python 0-indexed: 0,2,3,4,6,7
jj = np.array([1, 3, 4, 5, 7, 8]) - 1

T_list = []

# 3. First Figure: Plotting lines
plt.figure()
for j in range(len(jj)):
    # Directly grab the array and squeeze out the leading singleton dimension
    T_curr = s5c[jj[j]]["T"].squeeze()

    # This will now correctly output: Size of s5c(X).T: (29649,)
    print(f"Size of s5c({jj[j]+1}).T:", T_curr.shape)

    plt.plot(T_curr)
    T_list.append(T_curr)

plt.title("Line Plots")
plt.show()

# Convert the collected list into a 2D NumPy array
T = np.array(T_list)
print("Final T shape for pcolor:", T.shape)  # Should be (6, 29649)

# 4. Second Figure: Pseudocolor plot
plt.figure()
plt.pcolor(T, shading="flat")
plt.colorbar()
plt.title("pcolor(T)")
plt.show()

# 5. Print sizes of UV for all elements in s5c
print("\n--- UV Sizes ---")
for j in range(len(s5c)):
    UV_curr = s5c[j]["UV"].squeeze()
    print(f"Size of s5c({j+1}).UV:", UV_curr.shape)