#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 1

hloc = ["122"]
hloc_val = ["PC"]

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

t = 0
temp = ["rd"]
#suf = ["", "_L45", "_L65", "_L75", "_L85", "_L95"]
lat = ["0", "45", "65", "75", "85", "90"]
sgr = ["N", "A" , "D" , "C" , "B"]
#sgr = ["N", "R", "A" , "D" , "C" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 10, 25 , 50 , 100])*sgr_factor
#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor


melt_total = np.zeros((len(lat),len(sgr)))
h = 0
for v in range(len(lat)):
    for s in range(len(sgr)):
        if lat[v] == '75':
            fdir = f'{p_base}/{temp[t]}/{temp[t]}_{hloc[h]}{sgr[s]}'
        elif lat[v] == '0':
            fdir = f'{p_base}/rdn/rdn_{hloc[h]}{sgr[s]}'
        else:
            fdir = f'{p_base}/{temp[t]}L{lat[v]}/{temp[t]}/{temp[t]}_{hloc[h]}{sgr[s]}_L{lat[v]}'

        ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
        ds.load()
        melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
        #melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFlux.data)
        melt = melt / rho_fw * secPerYear

        dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
        dsMesh.load()
        areaCell = np.squeeze(dsMesh.areaCell.data)
        FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
        iii = FloatingMask == 1
        #melt_total[v, s] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
        melt_total[v, s] = np.sum(melt[iii] * areaCell[iii])
    melt_total[v, :] = melt_total[v, :] - melt_total[v, 0]

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit_lg = np.linspace(sgr_val[1], sgr_val[-1], 100)
xfit = np.linspace(0,sgr_val[-1],100)

plt.figure(figsize=(4, 4))
clr = 'rkbcmg'
smb = 'ops^h>'
#smb = 'ooooooo'

for v in range(len(lat)):
    plt.plot(sgr_val, melt_total[v,:], f'{clr[v]}:', linewidth=1, marker=f'{smb[v]}', fillstyle='full', markersize=4, label = f'Lat = {lat[v]}$^\circ$')
    popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[v, :])
    plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[v]}-', linewidth=1)
    popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[v, :], p0=10**10)
    plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[v]}--', linewidth=1)
#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
#plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.ylabel('Melt-flux anomaly (m$^3$/a)')
#if temp[t] == "rd":
#    plt.title('$T_b=1^\circ$C, $f=-1.409 \cdot 10^{-4}$ s$^{-1}$', fontsize = 8)
#else:
#    plt.title('$T_b=1^\circ$C, $f=0$ s$^{-1}$', fontsize=8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])
plt.title('$T_b=1^\circ$C', fontsize=8)

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_latitude.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

plt.figure(figsize=(4, 4))
for v in range(len(lat)):
    plt.loglog(sgr_val, melt_total[v,:], f'{clr[v]}:', linewidth=1, marker=f'{smb[v]}', fillstyle='full', markersize=4, label = f'Lat = {lat[v]}$^\circ$')
    popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[v, :])
    plt.loglog(xfit_lg, f_n_pow_2_3(xfit_lg, *popt), f'{clr[v]}-', linewidth=1)
    popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[v, :], p0=10**10)
    plt.loglog(xfit_lg, f_n_pow_1_3(xfit_lg, *popt), f'{clr[v]}--', linewidth=1)

plt.xlabel('$F_{s}$ (m$^3$/s)')
#plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.ylabel('Melt-flux anomaly (m$^3$/a)')
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
plt.title('$T_b=1^\circ$C', fontsize=8)

if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_latitude_lglg.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

