#!/usr/bin/env python3
import numpy as np
import xarray as xr
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import matplotlib.ticker as ticker
from scipy import signal


fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

'''
fname = 'lifw_tseries'
p_file = f'{fpath}{fname}.nc'

ds = xr.open_dataset(p_file)

# Step 4: Extract the variable '__xarray_dataarray_variable__'
variable_name = '__xarray_dataarray_variable__'
var = ds[variable_name]
print(var.shape)

time = np.arange(0,len(var))/12 + 2015
print(time.shape)

# Step 5: Plot the data (assuming variable is 1D and aligned with time)
plt.figure(figsize=(10, 6))
plt.plot(time, var, label=variable_name)
plt.xlabel('Time')
plt.ylabel(variable_name)
plt.title(f'{variable_name} over Time')
plt.grid(True)
plt.legend()
plt.show()
'''

fname = 'proxy_tseries'
fname_fismf = 'proxy_tseries_fismf'

ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Step 4: Extract the variable '__xarray_dataarray_variable__'
lifw = ds['lifw']
#ut = ds['ut']
utsq = ds['utsq']
utstar = ds['utstar']
time = np.arange(0,len(lifw))/12 + 2015

utsq_fismf = ds_fismf['utsq']

k = np.mean(lifw/utsq)

# Step 5: Plot the data (assuming variable is 1D and aligned with time)
plt.figure(figsize=(10, 6))
plt.plot(time, lifw, label='lifw')
#plt.plot(time, ut, label='ut')
plt.plot(time, utsq, label='utsq')
plt.plot(time, utstar, label='utstar')
#plt.plot(time, utsq_fismf, label='utsq_fismf')

plt.xlabel('Time')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
#plt.plot(time, ut/lifw, label='ut/lifw')
plt.plot(time, utsq/lifw, label='utsq/lifw')
plt.plot(time, utstar/lifw, label='utstar/lifw')
plt.xlabel('Time')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(time, lifw, label='lifw')
plt.plot(time, utsq/27, label='utsq/27')
plt.plot(time, utstar/27, label='utstar/27')
plt.xlabel('Time')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(time, lifw, label='lifw')
plt.plot(time, utsq*k, '--',label='utsq*k')
plt.plot(time, utsq_fismf*k, '--', label='utsq_fismf*k')
plt.xlabel('Time')
plt.grid(True)
plt.legend()
plt.show()