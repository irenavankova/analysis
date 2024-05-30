#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 1

h = 3
hloc = ["112", "122", "132" , "142"]
hloc_val = ["L", "PC", "PW" , "PE"]

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = ["rd", "rdn"]
sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
#sgr = ["N", "R", "A" , "D" , "C" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 25 , 50 , 75 , 100])*sgr_factor
#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor

melt_total = np.zeros((len(temp),len(sgr)))
for t in range(len(temp)):
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
        melt_total[t, s] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
    melt_total[t, :] = melt_total[t, :] - melt_total[t, 0]

#coef = numpy.polyfit(temp_val, melt_total, 2)
#tfit = np.linspace(-2, 4, 50)
#qfit = np.polyval(coef, tfit)

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit = np.linspace(0,sgr_val[-1],100)

plt.figure(figsize=(4, 4))
clr = 'rb'
smb = 's^'
rot = ["f", "f0"]
for t in range(len(temp)):
    plt.plot(sgr_val, melt_total[t,:], f'{clr[t]}{smb[t]}', fillstyle='none', label = f'{rot[t]} sim')
    #fit
    popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[t, :])
    plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1, label = '$a \cdot n^{1/3}$')
    popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[t, :])
    plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1, label = '$a \cdot n^{2/3}$')

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.title(f'$F_s$ location: {hloc_val[h]}')
plt.legend(loc=2, prop={'size': 6})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_fvsnf_{hloc_val[h]}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()


