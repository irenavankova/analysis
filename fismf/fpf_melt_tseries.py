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

opt_save = 1
fname = 'pANS_f701'
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

fname_xarr = 'lifw_reg_ave_tseries'
fname_pismf = 'utsq_reg_pismf_ave_tseries'
fname_fismf = 'utsq_reg_fismf_701_tseries'

#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_xarr = xr.open_dataset(f'{fpath}{fname_xarr}.nc')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Extract the DataArray (assuming variable is 'lifw')
#lifw = ds['lifw']  # dims: (Time, region)
lifwx = ds_xarr['lifw']  # dims: (Time, region)
n_time = len(ds_xarr['Time'])
time = np.arange(0,n_time)/12 + 2015

fs = 1/(time[1]-time[0])
fc = 1/((time[1]-time[0])*36)
def butter_lowpass_filter(data, fcc, fss, order=4):
    nyq = 0.5 * fss  # Nyquist frequency
    normal_cutoff = fcc / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f = ds_fismf['utsq']  # dims: (Time, region)

# Plot time series for each region
for region in lifwx.region.values:
    if region == "Antarctica":
        lifw_reg = lifwx.sel(region=region)
        utp_reg = utsq_p.sel(region=region)
        cAnt = np.mean(lifw_reg.values / utp_reg.values)
        print(cAnt)

# PLOT
#clr = ["black", "brown", "royalblue", "darkorange","lightskyblue"]
#FOR FILT:
clr = ["black", "darkorange","lightskyblue", "brown", "royalblue"]

#clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
#clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

#Lwide = np.array([1.25, 0.75, 0.75, 0.75, 0.75])
#smb = ["-", "-", "-", ":", ":"]
Lwide = np.array([0.5, 0.5, 0.5, 1.5, 1.5])
smb = ["-", "--", "--", "-", "-"]

iceshelves = ["Antarctica", "Bellingshausen", "Amundsen", "Ross", "East Antarctica", "Amery", "Dronning Maud Land", "Filchner-Ronne", "Larsens"]

ncols = 3
nrows = int(np.ceil((len(ds_xarr['region']))) / ncols)

fHeight = 10
fWidth = 18
cm = 1/2.54

fig, axes = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
#plt.subplots_adjust(hspace=cm*1.25)  # Increase vertical spacing
#plt.subplots_adjust(wspace=0.3)  # Increase vertical spacing

j = int(0); k = int(0)
for r, region in enumerate(lifwx.region.values):
    lifw_reg = lifwx.sel(region=region)
    utp_reg = utsq_p.sel(region=region)
    utf_reg = utsq_f.sel(region=region)
    c = np.mean(lifw_reg.values / utp_reg.values)

    utp_filt = butter_lowpass_filter(utp_reg, fc, fs)
    utf_filt = butter_lowpass_filter(utf_reg, fc, fs)

    print(c)
    s = 0
    axes[j, k].plot(time, lifw_reg, smb[s], color=clr[s], linewidth=Lwide[s]); s = s +1
    #axes[j, k].plot(time, utp_reg * c, smb[s], color=clr[s], linewidth=Lwide[s], label='utp_reg*k'); s = s +1
    #axes[j, k].plot(time, utf_reg * c, smb[s], color=clr[s], linewidth=Lwide[s], label='utf_reg*k'); s = s +1
    axes[j, k].plot(time, utp_reg * cAnt, smb[s], color=clr[s], linewidth=Lwide[s]); s = s +1
    axes[j, k].plot(time, utf_reg * cAnt, smb[s], color=clr[s], linewidth=Lwide[s]); s = s +1
    axes[j, k].plot(time, utp_filt * cAnt, smb[s], color=clr[s], linewidth=Lwide[s], label='EVOMELT'); s = s +1
    axes[j, k].plot(time, utf_filt * cAnt, smb[s], color=clr[s], linewidth=Lwide[s], label='FIXMELT')

    axes[j, k].set_title(region)
    axes[j, k].grid(True)

    fsize = 8
    axes[j, k].set_xlim(numpy.array([np.min(time), np.max(time)]))
    axes[j, k].set_title(iceshelves[r], fontsize=fsize - 1)
    axes[j, k].autoscale(enable=True, axis='both', tight=True)
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

    k = k + 1
    if k > ncols - 1:
        k = 0
        j = j + 1

fig.tight_layout()

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/Meltrates/Melt_tseries_{fname}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()




