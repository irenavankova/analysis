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

fname_xarr = 'lifw_reg_xarr_tseries'
fname_pismf = 'utsq_reg_ssp_tseries'
fname_fismf = 'utsq_reg_fismf_tseries'
fname_lifw_701 = 'lifw_reg_701_tseries'

#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_xarr = xr.open_dataset(f'{fpath}{fname_xarr}.nc')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')
ds_lifw_701 = xr.open_dataset(f'{fpath}{fname_lifw_701}.nc')


# Extract the DataArray (assuming variable is 'lifw')
#lifw = ds['lifw']  # dims: (Time, region)
lifwx = ds_xarr['lifw']  # dims: (Time, region)
n_time = len(ds_xarr['Time'])
time = np.arange(0,n_time)/12 + 2015
lifwx_701 = ds_lifw_701['lifw']  # dims: (Time, region)


utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f = ds_fismf['utsq']  # dims: (Time, region)


# Plot time series for each region
for region in lifwx.region.values:
    if region == "Antarctica":
        lifw_reg = lifwx.sel(region=region)
        utp_reg = utsq_p.sel(region=region)
        kAnt = np.mean(lifw_reg.values / utp_reg.values)
        print(kAnt)

for region in lifwx.region.values:
    plt.figure(figsize=(12, 6))

    lifw_reg = lifwx.sel(region=region)
    utp_reg = utsq_p.sel(region=region)
    utf_reg = utsq_f.sel(region=region)
    k = np.mean(lifw_reg.values / utp_reg.values)

    lifw_701 = lifwx_701.sel(region=region)

    print(k)
    plt.plot(time, lifw_reg, '-', label='lifw')
    plt.plot(time, utp_reg * k, '--', label='utp_reg*k')
    plt.plot(time, utf_reg * k, '--', label='utf_reg*k')
    plt.plot(time, utp_reg * kAnt, ':', label='utp_reg*kAnt')
    plt.plot(time, utf_reg * kAnt, ':', label='utf_reg*kAnt')
    plt.plot(time, lifw_701, '-', label='lifw_701')

    plt.title(region)
    plt.xlabel('Time')
    plt.grid(True)
    plt.legend()
    plt.show()
