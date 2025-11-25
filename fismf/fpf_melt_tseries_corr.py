#!/usr/bin/env python3
import numpy as np
import xarray as xr
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import matplotlib.ticker as ticker
from scipy.signal import butter, filtfilt
import iv_filt

opt_save = 1
opt_plot_raw = 1
opt_plot_ens = 0
ff = 0
#fname = 'pANS_f701'
fname = 'pANS_R2'
fpath = '/Users/irenavankova/Work/data_sim/FISMF/FISMF_E3SM_ouputs/ncfiles/post_derived/'

fname_lifw_hist = 'lifw_reg_hist_ave_tseries'
fname_lifw_pismf = 'lifw_reg_ave_tseries'
fname_lifw_fismf = 'lifw_reg_fismf_ave_tseries'

fname_hist = 'utsq_reg_hist_ave_tseries'
fname_pismf = 'utsq_reg_pismf_ave_tseries'
fname_fismf = 'utsq_reg_fismf_ave_tseries'

ensnum = ["701", "751","801"]


#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_lifw_hist = xr.open_dataset(f'{fpath}{fname_lifw_hist}.nc')
ds_lifw_pismf = xr.open_dataset(f'{fpath}{fname_lifw_pismf}.nc')
ds_lifw_fismf = xr.open_dataset(f'{fpath}{fname_lifw_fismf}.nc')

ds_hist = xr.open_dataset(f'{fpath}{fname_hist}.nc')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Extract the DataArray (assuming variable is 'lifw')
#lifw = ds['lifw']  # dims: (Time, region)
lifw_h = ds_lifw_hist['lifw']  # dims: (Time, region)
n_time = len(ds_lifw_hist['Time'])
time_h = np.arange(0,n_time)/12 + 1950

lifw_p = ds_lifw_pismf['lifw']  # dims: (Time, region)
n_time = len(ds_lifw_pismf['Time'])
time_ssp = np.arange(0,n_time)/12 + 2015

lifw_f = ds_lifw_fismf['lifw']  # dims: (Time, region)


utsq_h = ds_hist['utsq']  # dims: (Time, region)
utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f = ds_fismf['utsq']  # dims: (Time, region)

fs = 1/(time_ssp[1]-time_ssp[0])
fc = 1/((time_ssp[1]-time_ssp[0])*36)

# Plot time series for each region
for region in lifw_p.region.values:
    if region == "Antarctica":
        lifw_reg = lifw_p.sel(region=region)
        utp_reg = utsq_p.sel(region=region)
        kAnt = np.mean(lifw_reg.values / utp_reg.values)
        print(kAnt)

# PLOT
#clr = ["black", "brown", "royalblue", "darkorange","lightskyblue"]
#FOR FILT:
clr = ["black", "darkorange","lightskyblue", "brown", "royalblue"]

#clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
#clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

#Lwide = np.array([1.25, 0.75, 0.75, 0.75, 0.75])
#smb = ["-", "-", "-", ":", ":"]
#Lwide = np.array([0.5, 0.5, 0.5, 1.5, 1.5])
#smb = ["-", "--", "--", "-", "-"]
Lwide = 0.5
Lwide2 = 1.5

iceshelves = ["Antarctica", "Bellingshausen", "Amundsen", "Ross", "East Antarctica", "Amery", "Dronning Maud Land", "Filchner-Ronne", "Larsen"]
abc = 'abcdefghijklmnop'
ncols = 3
nrows = int(np.ceil((len(ds_lifw_pismf['region']))) / ncols)

fHeight = 14
fWidth = 16
cm = 1/2.54

fig, axes = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.4)  # Increase vertical spacing
plt.subplots_adjust(wspace=0.2)  # Increase horizontal spacing

j = int(0); k = int(0)
for r, region in enumerate(lifw_p.region.values):

    lh_reg = lifw_h.sel(region=region)
    lp_reg = lifw_p.sel(region=region)
    lf_reg = lifw_f.sel(region=region)

    uth_reg = utsq_h.sel(region=region)
    utp_reg = utsq_p.sel(region=region)
    utf_reg = utsq_f.sel(region=region)

    l_reg = np.append(lh_reg, lp_reg, axis=0)
    p_reg = np.append(uth_reg, utp_reg, axis=0)
    f_reg = np.append(uth_reg, utf_reg, axis=0)

    if opt_plot_raw == 1:
        #axes[j, k].plot(time_h, lh_reg, '-', color="black", linewidth=Lwide)
        axes[j, k].scatter(lh_reg, uth_reg * kAnt, marker=".", color="black", s=5)

        #axes[j, k].plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide)
        axes[j, k].scatter(lp_reg, utp_reg * kAnt, marker=".", color="black", s=5)

        #axes[j, k].plot(time_ssp, utf_reg * kAnt, '--', color="lightblue", linewidth=Lwide)



    axes[j, k].set_title(region)
    axes[j, k].grid(True)

    fsize = 8
    #axes[j, k].set_xlim(numpy.array([np.min(time_h), np.max(time_ssp)]))
    axes[j, k].set_title(f'{abc[r]}) {iceshelves[r]}', fontsize=fsize - 1)
    #axes[j, k].autoscale(enable=True, axis='both', tight=True)
    axes[j, k].tick_params(axis='both', labelsize=fsize)

    axes[j, k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')

    nloc = ticker.MaxNLocator(nbins=3)
    axes[j, k].yaxis.set_major_locator(nloc)
    tick_locations = nloc.tick_values(0, np.round(np.max(lp_reg)))
    tick_locations = np.round(tick_locations)

    # 4. Set the same tick locations for both X and Y axes
    axes[j, k].set_xticks(tick_locations)
    axes[j, k].set_yticks(tick_locations)

    if k == 0:
        axes[j, k].set_ylabel('$\dot{M}_d$ (Gt/a)', fontsize=fsize)

    if j == nrows-1:
        axes[j, k].set_xlabel('$\dot{M}_i$ (Gt/a)', fontsize=fsize)


    axes[j, k].set_aspect('equal', adjustable='datalim')  #
    axes[j, k].set_box_aspect(1)

    #if r == 0:
    #    axes[j, k].legend()
    #    axes[j, k].legend(fontsize=6)

    k = k + 1
    if k > ncols - 1:
        k = 0
        j = j + 1


#fig.tight_layout()

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/Meltrates/Melt_scatter_{fname}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()




