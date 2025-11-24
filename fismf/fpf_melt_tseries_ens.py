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
opt_plot_ens = 1
ff = 1
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

iceshelves = ["Antarctica", "Bellingshausen", "Amundsen", "Ross", "East Antarctica", "Amery", "Dronning Maud Land", "Filchner-Ronne", "Western Weddell"]
abc = 'abcdefghijklmnop'
ncols = 2
nrows = int(np.ceil((len(ds_lifw_pismf['region']))) / ncols)+1

fHeight = 17
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
        axes[j, k].plot(time_h, uth_reg * kAnt, '--', color="yellowgreen", linewidth=Lwide)

        #axes[j, k].plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide)
        axes[j, k].plot(time_ssp, utp_reg * kAnt, '--', color="navajowhite", linewidth=Lwide)

        axes[j, k].plot(time_ssp, utf_reg * kAnt, '--', color="lightblue", linewidth=Lwide)

    if opt_plot_ens == 1:
        for b in range(3):
            ds_pismf_ens = xr.open_dataset(f'{fpath}utsq_reg_pismf_{ensnum[b]}_tseries.nc')
            utsq_p_ens = ds_pismf_ens['utsq']  # dims: (Time, region)
            utp_ens = utsq_p_ens.sel(region=region)
            p_ens = np.append(uth_reg, utp_ens, axis=0)
            p_ens = iv_filt.butter_filter(p_ens, fc, fs, 'low')

            axes[j, k].plot(time_ssp, p_ens[len(time_h):] * kAnt, '-', color="darkorange", linewidth=Lwide+0.5)

            ds_fismf_ens = xr.open_dataset(f'{fpath}utsq_reg_fismf_{ensnum[b]}_tseries.nc')
            utsq_f_ens = ds_fismf_ens['utsq']  # dims: (Time, region)
            utf_ens = utsq_f_ens.sel(region=region)
            f_ens = np.append(uth_reg, utf_ens, axis=0)
            f_ens = iv_filt.butter_filter(f_ens, fc, fs, 'low')

            axes[j, k].plot(time_ssp, f_ens[len(time_h):] * kAnt, '-', color="cornflowerblue", linewidth=Lwide+0.5)



    if ff == 1:
        l_reg = iv_filt.butter_filter(l_reg, fc, fs, 'low')
        p_reg = iv_filt.butter_filter(p_reg, fc, fs, 'low')
        f_reg = iv_filt.butter_filter(f_reg, fc, fs, 'low')

        lh_reg = l_reg[0:len(time_h)]
        lp_reg = l_reg[len(time_h):]

        uth_reg = p_reg[0:len(time_h)]
        utp_reg = p_reg[len(time_h):]

        utf_reg = f_reg[len(time_h):]

        #axes[j, k].plot(time_h, lh_reg, '-', color="black", linewidth=Lwide2, label='lifw_p')
        axes[j, k].plot(time_h, uth_reg * kAnt, '-', color="darkolivegreen", linewidth=Lwide2, label='HIST')

        #axes[j, k].plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide2, label='lifw_p')
        axes[j, k].plot(time_ssp, utp_reg * kAnt, '-', color="brown", linewidth=Lwide2, label='EVOMELT')

        axes[j, k].plot(time_ssp, utf_reg * kAnt, '-', color="darkblue", linewidth=Lwide2, label='FIXMELT')

    axes[j, k].plot(time_ssp, lf_reg, '-', color="gold", linewidth=Lwide2, label='FIXMELT prescribed')

    axes[j, k].set_title(region)
    axes[j, k].grid(True)

    fsize = 8
    axes[j, k].set_xlim(numpy.array([np.min(time_h), np.max(time_ssp)]))
    axes[j, k].set_title(f'{abc[r]}) {iceshelves[r]}', fontsize=fsize - 1)
    #axes[j, k].autoscale(enable=True, axis='both', tight=True)
    axes[j, k].tick_params(axis='both', labelsize=fsize)

    axes[j, k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    axes[j, k].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))

    if k == 0:
        axes[j, k].set_ylabel('$\dot{M}$ (Gt/a)', fontsize=fsize)

    if j == nrows-1:
        axes[j, k].set_xlabel('Time (a)', fontsize=fsize)
    else:
        plt.setp(axes[j, k].get_xticklabels(), visible=False)


    if r == 0:
        axes[j, k].legend()
        axes[j, k].legend(fontsize=6)

    if r == 0:
        fig.delaxes(axes[0, 1])
        k = 0
        j = 1
    else:
        k = k + 1
        if k > ncols - 1:
            k = 0
            j = j + 1


#fig.tight_layout()

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/Meltrates/Melt_tseries_{fname}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()




