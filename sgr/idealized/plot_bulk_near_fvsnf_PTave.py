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

d_scale = 3000

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = ["rd", "rdn"]
sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
#sgr = ["N", "R", "A" , "D" , "C" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 25 , 50 , 75 , 100])*sgr_factor
#sgr_val = np.array([0, 1, 10, 25 , 50 , 100])*sgr_factor
hloc = [ "122", "132", "142"]
hloc_val = ["PC", "PW", "PE"]

fdir = f'{p_base}/{temp[0]}/{temp[0]}_{hloc[0]}R'
di = xarray.open_dataset(f'{fdir}/init.nc')
di.load()
x = np.squeeze(di.xCell.data)
y = np.squeeze(di.yCell.data)
z = np.squeeze(di.zCell.data)

for h in range(len(hloc)):
    fdir = f'{p_base}/{temp[0]}/{temp[0]}_{hloc[h]}R'
    df = xarray.open_dataset(f'{fdir}/sgr_data.nc')
    df.load()
    subglacialRunoffFlux = np.squeeze(df.subglacialRunoffFlux.data)
    sgr_iii = subglacialRunoffFlux > 0
    sgr_pt = np.where(sgr_iii)[0]

    dist2point = numpy.absolute(numpy.sqrt((x - x[sgr_pt]) ** 2 + (y - y[sgr_pt]) ** 2 + (z - z[sgr_pt]) ** 2))
    #ii = numpy.where(dist2point == dist2point.min())
    i_ave = dist2point < d_scale
    print(hloc_val[h])
    print(np.where(i_ave)[0])
    if hloc[h] == "122":
        i_ave_122 = i_ave
    elif hloc[h] == "132":
        i_ave_132 = i_ave
    elif hloc[h] == "142":
        i_ave_142 = i_ave

melt_total = np.zeros((len(temp), len(sgr), len(hloc)))
melt_ave = np.zeros((len(temp), len(sgr)))
for t in range(len(temp)):
    for h in range(len(hloc)):
        if hloc[h] == "122":
            i_av = i_ave_122
        elif hloc[h] == "132":
            i_av = i_ave_132
        elif hloc[h] == "142":
            i_av = i_ave_142
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
            i_float = FloatingMask == 1
            iii = np.logical_and(i_av, i_float)
            #print(hloc_val[h])
            #print(np.where(iii)[0])
            melt_total[t, s, h] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
        melt_total[t, :, h] = melt_total[t, :, h] - melt_total[t, 0, h]
    #print(np.shape(np.squeeze(np.sum(melt_total, axis=2))))
    #print(np.shape(melt_ave[t, :]))

melt_ave[:, :] = np.squeeze(np.sum(melt_total, axis=2)) / 3

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit = np.linspace(0,sgr_val[-1],100)

plt.figure(figsize=(4, 4))
clr = 'rb'
smb = 's^p'
rot = ["f", "f0"]

for t in range(len(temp)):
    #for h in range(len(hloc)):
        #plt.plot(sgr_val, np.squeeze(melt_total[t,:,h]), f'{clr[t]}{smb[h]}', fillstyle='none', label = f'{rot[t]} {hloc_val[h]}')
    plt.plot(sgr_val, melt_ave[t, :], f'{clr[t]}o', fillstyle='full', label=f'{rot[t]} ave')

    #fit
    popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_ave[t, :])
    plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[t]}--', linewidth=1, label = '$a \cdot n^{1/3}$')
    popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_ave[t, :])
    plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[t]}-', linewidth=1, label = '$a \cdot n^{2/3}$')


plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.title(f'$F_s$ location: Point')
plt.legend(loc=2, prop={'size': 6})
plt.grid()
plt.rcParams.update({'font.size': 8})

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_near_fvsnf_PTave.png', bbox_inches='tight', dpi=300)
else:
    plt.show()


