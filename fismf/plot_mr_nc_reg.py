#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt


ff = 1
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

fname_xarr = 'lifw_reg_ave_tseries'
fname_pismf = 'utsq_reg_pismf_ave_tseries'
fname_fismf_701 = 'utsq_reg_fismf_701_tseries'
fname_fismf_751 = 'utsq_reg_fismf_751_tseries'
fname_lifw_701 = 'lifw_reg_701_tseries'

#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_xarr = xr.open_dataset(f'{fpath}{fname_xarr}.nc')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
ds_fismf_701 = xr.open_dataset(f'{fpath}{fname_fismf_701}.nc')
ds_fismf_751 = xr.open_dataset(f'{fpath}{fname_fismf_751}.nc')
ds_lifw_701 = xr.open_dataset(f'{fpath}{fname_lifw_701}.nc')


# Extract the DataArray (assuming variable is 'lifw')
#lifw = ds['lifw']  # dims: (Time, region)
lifwx = ds_xarr['lifw']  # dims: (Time, region)
n_time = len(ds_xarr['Time'])
time = np.arange(0,n_time)/12 + 2015
lifwx_701 = ds_lifw_701['lifw']  # dims: (Time, region)

utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f_701 = ds_fismf_701['utsq']  # dims: (Time, region)
utsq_f_751 = ds_fismf_751['utsq']  # dims: (Time, region)

fs = 1/(time[1]-time[0])
fc = 1/((time[1]-time[0])*36)

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
    utf_701 = utsq_f_701.sel(region=region)
    utf_751 = utsq_f_751.sel(region=region)
    k = np.mean(lifw_reg.values / utp_reg.values)
    utf_ave = (utf_701 + utf_751)*0.5

    lifw_701 = lifwx_701.sel(region=region)

    # clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
    # clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

    print(k)
    Lwide = 1.25

    if ff == 1:
        lifw_reg = iv_filt.butter_filter(lifw_reg, fc, fs, 'low')
        utf_701 = iv_filt.butter_filter(utf_701, fc, fs, 'low')
        utf_751 = iv_filt.butter_filter(utf_751, fc, fs, 'low')
        utf_ave = iv_filt.butter_filter(utf_ave, fc, fs, 'low')
        lifw_701 = iv_filt.butter_filter(lifw_701, fc, fs, 'low')

    plt.plot(time, lifw_reg, '-', color="brown", linewidth=Lwide, label='lifw_ave')
    #plt.plot(time, utp_reg * k, '--', color="brown", linewidth=Lwide, label='utp_reg*k')
    #plt.plot(time, utf_701 * k, '--', color="royalblue", linewidth=Lwide, label='utf_701*k')
    #plt.plot(time, utf_751 * k, '--', color="darkolivegreen", linewidth=Lwide, label='utf_751*k')
    #plt.plot(time, utf_ave * k, '--', color="indigo", linewidth=Lwide, label='utf_ave*k')
    #plt.plot(time, utp_reg * kAnt, ':', color="darkorange", linewidth=Lwide, label='utp_reg*kAnt')
    plt.plot(time, utf_701 * kAnt, '--', color="lightskyblue", linewidth=Lwide, label='utf_701*kAnt')
    plt.plot(time, utf_751 * kAnt, '--', color="royalblue", linewidth=Lwide, label='utf_751*kAnt')
    plt.plot(time, utf_ave * kAnt, '-', color="indigo", linewidth=Lwide, label='utf_ave*kAnt')
    plt.plot(time, lifw_701, '--', color="darkorange", linewidth=Lwide, label='lifw_701')


    plt.title(region)
    plt.xlabel('Time')
    plt.grid(True)
    plt.legend()
    plt.show()
