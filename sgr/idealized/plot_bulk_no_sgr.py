#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = ["ra", "rb", "rc", "rd", "re", "rg"]
temp_val = np.array([-1.9, -1., 0, 1, 2 , 4])
sgr = 'R'
hloc = '112'

melt_total = np.zeros(len(temp))

for t in range(len(temp)):
    fdir = f'{p_base}/{temp[t]}/{temp[t]}_{hloc}{sgr}'

    ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
    ds.load()
    melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
    melt = melt/rho_fw*secPerYear

    dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
    dsMesh.load()
    areaCell = np.squeeze(dsMesh.areaCell.data)
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
    iii = FloatingMask == 1
    melt_total[t] = np.sum(FloatingMask[iii]*melt[iii]*areaCell[iii])/np.sum(areaCell[iii])


coef = numpy.polyfit(temp_val, melt_total, 2)
tfit = np.linspace(-2, 4, 50)
qfit = np.polyval(coef, tfit)
plt.figure(figsize=(4, 4))
plt.plot(temp_val, melt_total, 'ro', fillstyle = 'none' , label = 'sim')
plt.plot(tfit, qfit, 'k--' , linewidth=1, label = 'fit')
plt.xlabel('$T_b$ ($^\circ$C)')
plt.ylabel('Mean melt rate (m/a)')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()





