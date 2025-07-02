#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt



opt_save = 1
ff = 0
var_plot = 'siv'
if var_plot == 'sic':
    fac = 1e6
else:
    fac = 1e9

fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

ds = xr.open_dataset(f'{fpath}siv_max_over_cells.nc')
pmax = ds['pismf_max_over_cells']  # dims: (Time, region)
fmax = ds['fismf_max_over_cells']  # dims: (Time, region)
p10 = ds['pismf_10_over_cells']  # dims: (Time, region)
f10 = ds['fismf_10_over_cells']  # dims: (Time, region)
p2 = ds['pismf_2_over_cells']  # dims: (Time, region)
f2 = ds['fismf_2_over_cells']  # dims: (Time, region)
n_time = len(ds['Time'])
time = np.arange(0,n_time)/12 + 2015
plt.figure(figsize=(12, 6))
plt.plot(time, pmax, '-', color="darkorange", linewidth=1.5, label='pmax')
plt.plot(time, fmax, '-', color="lightskyblue", linewidth=1.5, label='fmax')
plt.plot(time, p2, '--', color="maroon", linewidth=1.5, label='p2')
plt.plot(time, f2, ':', color="indigo", linewidth=1.5, label='f2')
plt.plot(time, p10, '--', color="lightcoral", linewidth=1.5, label='p10')
plt.plot(time, f10, ':', color="dodgerblue", linewidth=1.5, label='f10')
plt.legend()
plt.show()

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
    plt.figure(figsize=(8, 4))

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

    fsize = 10
    plt.title(region,fontsize=fsize)
    plt.xlim(np.array([np.min(time), np.max(time)]))
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tick_params(axis='both', labelsize=fsize)
    plt.xlabel('Time (a)', fontsize=fsize)
    if var_plot == 'sic':
        plt.ylabel('Sea Ice Area (km^2)', fontsize=fsize)
    else:
        plt.ylabel('Sea Ice Volume (km^3)', fontsize=fsize)
    plt.grid(True)
    plt.legend()
    plt.legend(fontsize=6)
    if opt_save == 1:
        plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/SeaIceTseries/SIVtseries_{region.replace(" ", "_")}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()
