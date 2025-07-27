#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt

opt_save = 1
var_plot = 'tau_zonal'
#var_plot = 'tau_merid'

fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

fname_pismf = 'tau_reg_pismf_ave_tseries'
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
fname_fismf = 'tau_reg_fismf_701_tseries'
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Extract the DataArray
sip = ds_pismf[var_plot]  # dims: (Time, region)
sif = ds_fismf[var_plot]  # dims: (Time, region)

n_time = len(ds_pismf['Time'])
time = np.arange(0,n_time)/12 + 2015

fs = 1/(time[1]-time[0])
fc = 1/((time[1]-time[0])*12*3)

for region in sip.region.values:
    plt.figure(figsize=(8, 4))

    sip_reg = sip.sel(region=region)
    sif_reg = sif.sel(region=region)

    # clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
    # clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

    Lwide = 1.25

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
    if var_plot == 'tau_zonal':
        plt.ylabel('Zonal windstress (N/m^2)', fontsize=fsize)
    else:
        plt.ylabel('Meridional windstress (N/m^2)', fontsize=fsize)
    plt.grid(True)
    plt.legend()
    plt.legend(fontsize=6)
    if opt_save == 1:
        plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/TauTseries/{var_plot}_{region.replace(" ", "_")}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()
