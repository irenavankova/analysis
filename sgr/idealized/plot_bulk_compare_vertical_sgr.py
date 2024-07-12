#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 1

hloc = ["111", "112", "110"]
hloc_val = ["B", "U", "T"]
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_data_conserve_04_yesC'

temp = ["f100", "f101", "f102"]
sgr = ["N", "", "A", "C", "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 50, 100])*sgr_factor
temp_val = np.array([1, 0, -1])

#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor

melt_total = np.zeros((len(temp), len(hloc), len(sgr)))
for t in range(len(temp)):
    for h in range(len(hloc)):
        for s in range(len(sgr)):
            fdir = f'{p_base}/{temp[t]}/{temp[t]}_{hloc[h]}000{sgr[s]}'

            ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
            ds.load()
            #melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
            melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFlux.data)
            melt = melt / rho_fw * secPerYear

            dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
            dsMesh.load()
            areaCell = np.squeeze(dsMesh.areaCell.data)
            FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
            iii = FloatingMask == 1
            melt_total[t, h, s] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
        melt_total[t, h, :] = melt_total[t, h, :] - melt_total[t, h, 0]

plt.figure(figsize=(4, 4))
clr = 'rbk'
smb = 's^o'
ltp = ["-", "--", ":"]
for t in range(len(temp)):
    for h in range(len(hloc)):
        plt.plot(sgr_val, np.squeeze(melt_total[t,h,:]), f'{clr[t]}{ltp[h]}', linewidth=1, marker=f'{smb[h]}', fillstyle='none', markersize=4, label = f'{hloc_val[h]}, $T_b$={temp_val[t]}$^\circ$C')

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.title('$f=-1.409 \cdot 10^{-4}$ s$^{-1}$', fontsize = 8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_vertical_sgr.png', bbox_inches='tight', dpi=300)
else:
    plt.show()


