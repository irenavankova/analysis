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
from scipy.optimize import curve_fit


opt_save = 1

Tmax = 50
fyrs = f'y{Tmax}'

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/melt_flux_tseries_{fyrs}.nc'
dsOut = xarray.open_dataset(p_file)

MeltFlux = np.array(dsOut.MeltFlux.data)
Time = np.array(dsOut.time.data) + 1 - 1/24

print('Time')
print(Time[0])
print(Time[11])
print(Time[-1])

iceshelves = np.array(dsOut.iceshelves.data)
sims = np.array(dsOut.sims.data)

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/SgrMassFlux.nc'
dsOut = xarray.open_dataset(p_file)
SgrMassFlux = np.array(dsOut.SgrMassFlux.data)
IceShelfArea = np.array(dsOut.IceShelfArea.data)

#kgint = 1e3
#SgrMassFlux = SgrMassFlux/kgint
SgrMassFlux = SgrMassFlux/(IceShelfArea)
#sgr_unit = 'kg$\cdot$m$^{-2}\cdot$a$^{-1}$'
sgr_unit = 'kg/(m$^{2}$a)'

#ii = np.array([0, 12, 8, 7, 9, 17, 10, 1, 2, 5, 16, 15])
#ii = np.arange(0,23,1)
#ipl = ["Antarctica", "Antarctica-Amery", "George_VI", "Pine_Island", "Thwaites", "Getz", "NSS", "Ross", "Totten", "Shackleton", "Amery", "RoiB", "Fimbul", "Riiser-Larsen/Brunt", "Filchner-Ronne","Larsen_C"]
#ipl = ["Antarctica", "Antarctica-Amery", "George_VI","Thwaites", "Getz", "NSS", "Ross"]
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

Tmin = 20

itmin = Tmin*12
itmax = Tmax*12

#ts1 = np.array([20, 21, 21, 24, 31, 36, 41, 46])
#ts2 = np.array([21, 23, 26, 32, 40, 40, 50, 50])
#ts1 = np.array([20, 21, 21, 24, 31, 41])
#ts2 = np.array([21, 23, 26, 32, 40, 50])
#ts1 = np.array([21, 24, 31, 41])
#ts2 = np.array([23, 32, 40, 50])#

#ORIGINAL SUBMISSION:
#ts1 = np.array([21, 41])
#ts2 = np.array([23, 50])

#ts1 = np.array([21, 41]) #includes years 21, 22, 23
#ts2 = np.array([24, 51]) #includes years 41, 42, ..., 50
#tseg = f'y3'

#ts1 = np.array([22, 41]) #includes years 22, 23
#ts2 = np.array([24, 51]) #includes years 41, 42, ..., 50
#tseg = f'y2'

#ts1 = np.array([21, 41]) #includes years 21
#ts2 = np.array([22, 51]) #includes years 41, 42, ..., 50
#tseg = f'y1'

ts1 = np.array([22, 41]) #includes years 22, 22
ts2 = np.array([23, 51]) #includes years 41, 42, ..., 50
tseg = f'main'

#ts1 = np.array([20])
#ts2 = np.array([22])
#tseg = f'{ts1[0]}_{ts2[0]}'


#clr = ["brown", "orange", "deepskyblue" , "black"]
#clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
#clr = ["royalblue", "darkorange", "lightskyblue", "moccasin", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
clr = ["black", "red", "lightskyblue", "moccasin", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
abc = 'abcdefghijklmnopq'
ncols = 2
#nrows = (len(ii)+1) // ncols+1
nrows = (len(ii)+1) // ncols

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

def f_n_lin(x, a):
    return a * x

fHeight = 20
fWidth = 12
cm = 1/2.54

fig, axes = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=cm*1.25,wspace=0.75*cm)  # Increase vertical spacing


#np.shape(H)[1]
j = 0; k = 0
#for ind_r in range(np.shape(MeltFlux)[2]):
for ind_r in range(len(ii)):
    print('SGR mass flux')
    print(SgrMassFlux[ii[ind_r]])
    sgr_fit = np.arange(0, 9, 0.1)#*SgrMassFlux[ii[ind_r]]
    sgr = np.array([0, 1, 4, 8])#*SgrMassFlux[ii[ind_r]]

    for ind_ts in range(len(ts1)):
        melt = np.zeros(len(sims))
        ind_t = np.where((Time > ts1[ind_ts]) & (Time < ts2[ind_ts]))
        print('Time')
        print(ind_t[0][0])
        print(Time[ind_t[0][0]])
        print(ind_t[-1][0])
        print(Time[ind_t[-1][-1]])
        #print(ind_t[-1])

        for s in range(len(sims)):
            yplt = (np.squeeze(MeltFlux[:, s, ii[ind_r]]))

            melt[s] = np.mean(yplt[ind_t])

        axes[j,k].plot(sgr, melt-melt[0] ,'*', color = clr[ind_ts] ,linewidth=0.75, label = sims[s], markersize=6)

        popt, pcov = curve_fit(f_n_pow_1_3, sgr, melt-melt[0], p0=10 ** 10)
        axes[j,k].plot(sgr_fit, f_n_pow_1_3(sgr_fit, *popt), '--', color = clr[ind_ts] ,linewidth=0.75)
        popt, pcov = curve_fit(f_n_pow_2_3, sgr, melt-melt[0])
        axes[j,k].plot(sgr_fit, f_n_pow_2_3(sgr_fit, *popt), '-', color = clr[ind_ts] ,linewidth=0.75)
        popt, pcov = curve_fit(f_n_lin, sgr, melt - melt[0])
        axes[j, k].plot(sgr_fit, f_n_lin(sgr_fit, *popt), ':', color=clr[ind_ts], linewidth=0.75)

    fsize = 8
    axes[j,k].set_title(f'{abc[ind_r]}) {npl[ind_r]}, $F_s$/$A$ = {round(SgrMassFlux[ii[ind_r]])} {sgr_unit}', fontsize=fsize-1)
    axes[j,k].autoscale(enable=True, axis='both', tight=True)
    axes[j,k].tick_params(axis='both', labelsize=fsize)
    if k == 0:
        axes[j, k].set_ylabel('$\Delta \dot{M}$ (Gt/a)', fontsize=fsize)

    if j == nrows-1:
        axes[j, k].set_xlabel('$F_s$ factor', fontsize=fsize)
    else:
        plt.setp(axes[j, k].get_xticklabels(), visible=False)
    axes[j,k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    axes[j, k].yaxis.set_major_locator(ticker.MaxNLocator(nbins = 3))
    axes[j, k].set_xlim([-0.5, 8.5])
    ymin, ymax = axes[j, k].get_ylim()
    y_offset = 0.1 * max(abs(ymin), abs(ymax))
    axes[j, k].set_ylim(ymin-y_offset, ymax + y_offset)
    axes[j, k].set_xticks([0, 2, 4, 6, 8])

    k = k + 1
    if k > ncols-1:
        k = 0
        j = j + 1

#fig.delaxes(axes[nrows-1][ncols-2])
#fig.delaxes(axes[nrows-1][ncols-1])

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/melt_scaling/melt_12scaling_{tseg}_R2.png', bbox_inches='tight', dpi=300)
else:
    plt.show()