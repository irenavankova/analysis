#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt


ff = 0
var_plot = 'siv'
if var_plot == 'sic':
    fac = 1e6
else:
    fac = 1e9

fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

ds = xr.open_dataset(f'{fpath}{'siv_max_over_cells'}.nc')
sip = ds_pismf[var_plot]  # dims: (Time, region)


fname_pismf = 'si_reg_pismf_ave_tseries'
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
fname_fismf = 'si_reg_fismf_701_tseries'
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Extract the DataArray
sip = ds_pismf[var_plot]  # dims: (Time, region)
sif = ds_fismf[var_plot]  # dims: (Time, region)

n_time = len(ds_pismf['Time'])
time = np.arange(0,n_time)/12 + 2015

fs = 1/(time[1]-time[0])
fc = 1/((time[1]-time[0])*36)

for region in sip.region.values:
    plt.figure(figsize=(12, 6))

    sip_reg = sip.sel(region=region)/fac
    sif_reg = sif.sel(region=region)/fac

    # clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
    # clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

    Lwide = 1.25

    if ff == 1:
        lifw_reg = iv_filt.butter_filter(lifw_reg, fc, fs, 'low')
        utf_701 = iv_filt.butter_filter(utf_701, fc, fs, 'low')
        utf_751 = iv_filt.butter_filter(utf_751, fc, fs, 'low')
        utf_ave = iv_filt.butter_filter(utf_ave, fc, fs, 'low')
        lifw_701 = iv_filt.butter_filter(lifw_701, fc, fs, 'low')

    plt.plot(time, sip_reg, '-', color="darkorange", linewidth=Lwide, label='pismf')
    plt.plot(time, sif_reg, '--', color="lightskyblue", linewidth=Lwide, label='fismf_701')
    sip_reg = iv_filt.butter_filter(sip_reg, fc, fs, 'low')
    sif_reg = iv_filt.butter_filter(sif_reg, fc, fs, 'low')
    plt.plot(time, sip_reg, '-', color="brown", linewidth=2.5, label='pismf')
    plt.plot(time, sif_reg, '-', color="royalblue", linewidth=2.5, label='fismf_701')

    #plt.plot(time, utp_reg * k, '--', color="brown", linewidth=Lwide, label='utp_reg*k')
    #plt.plot(time, utf_701 * k, '--', color="royalblue", linewidth=Lwide, label='utf_701*k')
    #plt.plot(time, utf_751 * k, '--', color="darkolivegreen", linewidth=Lwide, label='utf_751*k')
    #plt.plot(time, utf_ave * k, '--', color="indigo", linewidth=Lwide, label='utf_ave*k')
    #plt.plot(time, utp_reg * kAnt, ':', color="darkorange", linewidth=Lwide, label='utp_reg*kAnt')
    #plt.plot(time, utf_701 * kAnt, '--', color="lightskyblue", linewidth=Lwide, label='utf_701*kAnt')
    #plt.plot(time, utf_751 * kAnt, '--', color="royalblue", linewidth=Lwide, label='utf_751*kAnt')
    #plt.plot(time, utf_ave * kAnt, '-', color="indigo", linewidth=Lwide, label='utf_ave*kAnt')
    #plt.plot(time, lifw_701, '--', color="darkorange", linewidth=Lwide, label='lifw_701')

    plt.title(region)
    plt.xlabel('Time')
    if var_plot == 'sic':
        plt.ylabel('Sea Ice Area (km^2)')
    else:
        plt.ylabel('Sea Ice Volume (km^3)')
    plt.grid(True)
    plt.legend()
    plt.show()
