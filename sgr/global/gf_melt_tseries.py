#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import matplotlib.ticker as ticker
from scipy import signal

opt_plot = 0
opt_save = 1

if opt_plot == 0:
    fname = 'val'
elif opt_plot == 1:
    fname = 'dif'
elif opt_plot == 2:
    fname = 'perc'
elif opt_plot == 3:
    fname = 'var'

Tmax = 50
fyrs = f'y{Tmax}'

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/melt_flux_tseries_{fyrs}.nc'
dsOut = xarray.open_dataset(p_file)
MeltFlux = np.array(dsOut.MeltFlux.data)
Time = np.array(dsOut.time.data)
iceshelves = np.array(dsOut.iceshelves.data)
sims = np.array(dsOut.sims.data)

Tmax2 = 110
fyrs = f'y{Tmax2}'

p_file2 = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/melt_flux_tseries_{fyrs}.nc'
dsOut2 = xarray.open_dataset(p_file2)
MeltFlux2 = np.array(dsOut2.MeltFlux.data)
Time2 = np.array(dsOut2.time.data)

#ii = np.array([0, 19, 12, 8, 7, 9, 17, 10, 1, 2, 5, 18, 16, 15])
#ii = np.arange(0,23,1)
#ipl = ["Antarctica", "Antarctica-Amery", "George_VI", "Pine_Island", "Thwaites", "Getz", "NSS", "Ross", "Totten", "Shackleton", "Amery", "RoiB", "Fimbul", "Riiser-Larsen/Brunt", "Filchner-Ronne","Larsen_C"]
ipl = ["Antarctica", "Antarctica-Amery", "Pine_Island", "Thwaites", "Getz", "NSS", "Ross", "Totten", "Amery", "Fimbul", "Filchner-Ronne","Larsen_C"]
npl = ["Antarctica", "Antarctica-Amery", "Pine Island", "Thwaites", "Getz", "NSS", "Ross", "Totten", "Amery", "Fimbul", "Filchner-Ronne","Larsen C"]

indices = []
for i, item1 in enumerate(iceshelves):
    if item1 in ipl:
        indices.append(i)

ii = indices
print(ii)

print(iceshelves)
print(iceshelves[ii])

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/SgrMassFlux.nc'
dsOut = xarray.open_dataset(p_file)
SgrMassFlux = np.array(dsOut.SgrMassFlux.data)

SgrMassFlux = SgrMassFlux/1e12
sgr_unit = 'GT/a'

Tmin = 20
Toffset = 20
itmin = Tmin*12
itmax = Tmax*12
Time = Time-Toffset

itmax2 = Tmax2*12
Time2 = Time2-Toffset


#clr = ["brown", "orange", "deepskyblue" , "black"]
clr = ["brown", "darkorange", "deepskyblue", "indigo"]


ncols = 2
nrows = (len(ii)+1) // ncols

#fHeight = 20
fHeight = 16
fWidth = 16
cm = 1/2.54

fig, axes = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=cm*1.25)  # Increase vertical spacing


#np.shape(H)[1]
j = 0; k = 0
#for ind_r in range(np.shape(MeltFlux)[2]):
for ind_r in range(len(ii)):
    print(SgrMassFlux[ii[ind_r]])

    for s in range(len(sims)):

        if s == 0:
            #smb = 'k:'
            smb = '-'
            if opt_plot > 0:
                ref = np.squeeze(MeltFlux[itmin:itmax, s, ii[ind_r]])
            else:
                ref = np.squeeze(MeltFlux2[itmin:itmax2, s, ii[ind_r]])
        else:
            smb = '-'

        #print(np.squeeze(MeltFlux[:, s, ind_r]))
        if (s > 1) | (opt_plot > 0):
            yplt = (np.squeeze(MeltFlux[itmin:itmax, s, ii[ind_r]]))
            tplt = Time[itmin:itmax]
        else:
            yplt = (np.squeeze(MeltFlux2[itmin:itmax2, s, ii[ind_r]]))
            tplt = Time2[itmin:itmax2]

        if opt_plot == 0:
            axes[j,k].plot(tplt, yplt ,smb, color = clr[s] ,linewidth=0.75, label = sims[s])
        elif opt_plot == 1:
            axes[j,k].plot(tplt, yplt-ref ,smb, color = clr[s] ,linewidth=0.75, label = sims[s])
        elif opt_plot == 2:
            axes[j,k].plot(tplt, yplt/ref ,smb, color = clr[s] ,linewidth=0.75, label = sims[s])
        elif opt_plot == 3:
            axes[j,k].plot(tplt, signal.detrend(yplt) ,smb, color = clr[s] ,linewidth=0.75, label = sims[s])

    fsize = 8
    axes[j,k].set_xlim(numpy.array([Tmin, Tmax2])-Toffset)
    #axes[j,k].set_xlabel('Time (a)', fontsize=fsize)
    #axes[j,k].set_ylabel('Integrated melf flux (GT/a)', fontsize=fsize)

    #axes[j,k].set_legend(loc = 2,fontsize=12)
    #axes[j,k].set_title(iceshelves[ii[ind_r]], fontsize=fsize)
    axes[j,k].set_title(f'{npl[ind_r]}, $F_s$ = {round(SgrMassFlux[ii[ind_r]],1)} {sgr_unit}', fontsize=fsize-1)
    #axes[j,k].set_title(f'{iceshelves[ii[ind_r]]}, $F_s$ = {round(SgrMassFlux[ii[ind_r]],1)} {sgr_unit}', fontsize=fsize-1)
    axes[j,k].autoscale(enable=True, axis='both', tight=True)
    axes[j,k].tick_params(axis='both', labelsize=fsize)
    if k == 0:
        if opt_plot == 0:
            axes[j, k].set_ylabel('$\dot{M}$ (GT/a)', fontsize=fsize)
        elif opt_plot == 1:
            axes[j, k].set_ylabel('$\dot{M}$ diff (GT/a)', fontsize=fsize)
        elif opt_plot == 2:
            axes[j, k].set_ylabel('$\dot{M} perc$', fontsize=fsize)
        elif opt_plot == 3:
            axes[j, k].set_ylabel('$\delta \dot{M}$ (GT/a)', fontsize=fsize)
    if j == nrows-1:
        axes[j, k].set_xlabel('Time (a)', fontsize=fsize)
    else:
        plt.setp(axes[j, k].get_xticklabels(), visible=False)
    axes[j,k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    axes[j, k].yaxis.set_major_locator(ticker.MaxNLocator(nbins = 3))

    if npl[ind_r] == 'Filchner-Ronne':
        ymin, ymax = axes[j,k].get_ylim()
        axes[j, k].set_ylim(numpy.array([-0.1, ymax]))

    k = k + 1
    if k > ncols-1:
        k = 0
        j = j + 1

if opt_save == 1:
    #plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/melt_tseries/melt_12tseries_{fyrs}_{fname}.png', bbox_inches='tight', dpi=300)
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/melt_tseries/melt_12tseries_yall_{fname}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()