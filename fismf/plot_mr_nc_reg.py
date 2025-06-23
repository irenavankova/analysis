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

fname = 'lifw_reg_tseries'
fname_xarr = 'lifw_reg_xarr_tseries'

ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_xarr = xr.open_dataset(f'{fpath}{fname_xarr}.nc')

# Extract the DataArray (assuming variable is 'lifw')
lifw = ds['lifw']  # dims: (Time, region)
lifwx = ds_xarr['lifw']  # dims: (Time, region)
n_time = len(ds['Time'])
time = np.arange(0,n_time)/12 + 2015


# Plot time series for each region
for region in lifw.region.values:
    plt.figure(figsize=(12, 6))
    # Select the time series for the region
    ts = lifw.sel(region=region)
    plt.plot(time, ts, label=region)
    tsx = lifwx.sel(region=region)
    plt.plot(time, tsx, '--', label=region)
    plt.xlabel('Time')
    plt.grid(True)
    plt.legend()
    plt.show()
