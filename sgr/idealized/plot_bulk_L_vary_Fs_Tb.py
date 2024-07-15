#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 0

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = ["ra", "rb", "rc", "rd", "re", "rg"]
temp_val = np.array([-1.9, -1, 0, 1, 2, 4])
sgr = ["N", "R", "A", "C", "B"]
#sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
sgr_factor = 7.2169
sgr_factor_sd = np.round(sgr_factor, 1)
sgr_val = np.array([0, 1, 10 , 50 , 100])*sgr_factor
sgr_val_sd = sgr_val/sgr_factor*sgr_factor_sd
#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor

melt_total = np.zeros((len(temp),len(sgr)))
meltTot_total = np.zeros((len(temp),len(sgr)))
for t in range(len(temp)):
    for s in range(len(sgr)):
        fdir = f'{p_base}/{temp[t]}/{temp[t]}_112{sgr[s]}'

        ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
        ds.load()
        meltTot = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
        meltTot = meltTot / rho_fw * secPerYear
        melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFlux.data)
        melt = melt / rho_fw * secPerYear

        dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
        dsMesh.load()
        areaCell = np.squeeze(dsMesh.areaCell.data)
        FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
        iii = FloatingMask == 1
        melt_total[t, s] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
        meltTot_total[t, s] = np.sum(FloatingMask[iii] * meltTot[iii] * areaCell[iii]) / np.sum(areaCell[iii])
    melt_total[t, :] = melt_total[t, :] - melt_total[t, 0]
    meltTot_total[t, :] = meltTot_total[t, :] - meltTot_total[t, 0]

#coef = numpy.polyfit(temp_val, melt_total, 2)
#tfit = np.linspace(-2, 4, 50)
#qfit = np.polyval(coef, tfit)

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit = np.linspace(0,sgr_val[-1],100)

plt.figure(figsize=(4, 4))
clr = 'kbrcgy'
smb = 'o^spdh'
for t in range(len(temp)):
    plt.plot(sgr_val, melt_total[t,:], f'{clr[t]}:', linewidth=1, marker='o', fillstyle='none', markersize=4)
    plt.plot(sgr_val, meltTot_total[t,:], f'{clr[t]}', linewidth=1, linestyle='none', marker='x', fillstyle='none', markersize=4)
    plt.plot(sgr_val, meltTot_total[t,:], f'{clr[t]}--', linewidth=1, label = f'$T_b$ = {temp_val[t]}$^\circ$C')

    #fit
    #popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[t, :])
    #plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1, label = '$a \cdot n^{1/3}$')
    #popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[t, :])
    #plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1, label = '$a \cdot n^{2/3}$')

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.legend(loc=2, prop={'size': 6})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_L_vary_Fs_Tb_Fax.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

plt.figure(figsize=(4, 4))
clr = 'kbrcgy'
smb = 'o^spdh'
for s in range(1, len(sgr)):
    plt.plot(temp_val, melt_total[:,s], f'{clr[s]}:', linewidth=1, marker='o', fillstyle='none', markersize=4)
    plt.plot(temp_val, meltTot_total[:,s], f'{clr[s]}', linewidth=1, linestyle='none', marker='x', fillstyle='none', markersize=4)
    plt.plot(temp_val, meltTot_total[:,s], f'{clr[s]}--', linewidth=1, label = f'$F_s$ = {sgr_val_sd[s]}m$^3$/s')

    #fit
    #popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[t, :])
    #plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1, label = '$a \cdot n^{1/3}$')
    #popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[t, :])
    #plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1, label = '$a \cdot n^{2/3}$')

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$T_b$ ($^\circ$C)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.legend(loc=2, prop={'size': 6})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_L_vary_Fs_Tb_Tax.png', bbox_inches='tight', dpi=300)
else:
    plt.show()


