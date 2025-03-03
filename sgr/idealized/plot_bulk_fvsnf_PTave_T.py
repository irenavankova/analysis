#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 0
lglg = 0

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = ["rd", "rdn"]
sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
#sgr = ["N", "R", "A" , "D" , "C" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 25 , 50 , 75 , 100])*sgr_factor
#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor
hloc = [ "122", "132", "142"]
hloc_val = ["PC", "PW", "PE"]

melt_total = np.zeros((len(temp), len(sgr), len(hloc)))
melt_ave = np.zeros((len(temp), len(sgr)))

tmax_total = np.zeros((len(temp), len(sgr), len(hloc)))
tmax_ave = np.zeros((len(temp), len(sgr)))
for t in range(len(temp)):
    for h in range(len(hloc)):
        for s in range(len(sgr)):
            fdir = f'{p_base}/{temp[t]}/{temp[t]}_{hloc[h]}{sgr[s]}'

            ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
            ds.load()
            melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
            melt = melt / rho_fw * secPerYear

            dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
            dsMesh.load()
            areaCell = np.squeeze(dsMesh.areaCell.data)
            FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
            iii = FloatingMask == 1
            melt_total[t, s, h] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])

            maxLevelCell = dsMesh.maxLevelCell.data - 1
            temp_clim = np.squeeze(ds.timeMonthly_avg_activeTracers_temperature.data)
            ix = np.where(iii)[0]
            temp_clim = temp_clim[ix, maxLevelCell[ix]]
            #tmax_total[t, s, h] = np.mean(temp_clim, axis=0)
            tmax_total[t, s, h] = np.max(temp_clim, axis=0)
        melt_total[t, :, h] = melt_total[t, :, h] - melt_total[t, 0, h]
        #tmax_total[t, :, h] = tmax_total[t, :, h] - tmax_total[t, 0, h]
    #print(np.shape(np.squeeze(np.sum(melt_total, axis=2))))
    #print(np.shape(melt_ave[t, :]))

melt_ave[:, :] = np.squeeze(np.sum(melt_total, axis=2)) / 3
tmax_ave[:, :] = np.squeeze(np.sum(tmax_total, axis=2)) / 3

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

if lglg == 1:
    xfit = np.linspace(sgr_val[1], sgr_val[-1], 100)
else:
    xfit = np.linspace(0,sgr_val[-1],100)

plt.figure(figsize=(4, 4))
clr = 'rb'
smb = 's^p'
rot = ["rotating", "non-rotating"]

for t in range(len(temp)):
    #for h in range(len(hloc)):
        #plt.plot(sgr_val, np.squeeze(melt_total[t,:,h]), f'{clr[t]}{smb[h]}', fillstyle='none', label = f'{rot[t]} {hloc_val[h]}')
    if lglg == 1:
        plt.loglog(sgr_val, melt_ave[t, :], f'{clr[t]}o', fillstyle='full', markersize=4, label=f'{rot[t]}')
    else:
        plt.plot(sgr_val, melt_ave[t, :], f'{clr[t]}o', fillstyle='full', markersize=4, label=f'{rot[t]}')

    #fit
    popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_ave[t, :])
    if lglg == 1:
        plt.loglog(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1)
    else:
        plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1)
    popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_ave[t, :])
    if lglg == 1:
        plt.loglog(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1)
    else:
        plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1)


plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.title(f'Channelized, all')
plt.legend(loc=2, prop={'size': 8})
if lglg == 0:
    plt.grid()
plt.rcParams.update({'font.size': 8})

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk/fvsnf'
if opt_save == 1:
    if lglg == 1:
        plt.savefig(f'{dir_fig_save}/plot_bulk_fvsnf_PTave_loglog.png', bbox_inches='tight', dpi=300)
    else:
        plt.savefig(f'{dir_fig_save}/plot_bulk_fvsnf_PTave.png', bbox_inches='tight', dpi=300)
else:
    plt.show()


plt.figure(figsize=(4, 4))

smb = 'o^sphdx'
for t in range(len(temp)):
    for h in range(len(hloc)):
        plt.plot(sgr_val, np.squeeze(tmax_total[t, :, h]), f'{clr[t]}{smb[h]}', fillstyle='none', markersize=4, label=f'{rot[t]}')
        #plt.plot(sgr_val, np.squeeze(tmax_total[t, :, h]))
plt.show()

