#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 1

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
        melt_total[t, s] = np.sum(melt[iii] * areaCell[iii])
        meltTot_total[t, s] = np.sum(meltTot[iii] * areaCell[iii])
        #melt_total[t, s] = np.sum(FloatingMask[iii] * melt[iii] * areaCell[iii]) / np.sum(areaCell[iii])
        #meltTot_total[t, s] = np.sum(FloatingMask[iii] * meltTot[iii] * areaCell[iii]) / np.sum(areaCell[iii])
    melt_total[t, :] = melt_total[t, :] - melt_total[t, 0]
    meltTot_total[t, :] = meltTot_total[t, :] - meltTot_total[t, 0]

#melt_fit = melt_total
melt_fit = meltTot_total

print(sgr_val.shape)
print(temp_val.shape)
print(melt_fit.shape)

tM, sM = np.meshgrid(temp_val, sgr_val, indexing='ij')

print(sM.shape)
print(tM.shape)
print(melt_total.shape)

t = np.squeeze(tM.reshape((len(sgr_val)*len(temp_val)), 1))
s = np.squeeze(sM.reshape((len(sgr_val)*len(temp_val)), 1))
m = np.squeeze(melt_fit.reshape((len(sgr_val)*len(temp_val)), 1))

print(t.shape)

def func(ts, k, To):
    t, s = ts
    return k*(t - To) * s**(2/3)
    #return k * (t - To) + s


# Perform curve fitting
popt, pcov = curve_fit(func, (t, s), m)

k_fit = popt[0]
To_fit = popt[1]

print(popt)

t_range = np.linspace(temp_val[0], temp_val[-1], 50)
s_range = np.linspace(sgr_val[0], sgr_val[-1], 50)
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(sgr_val[1:len(sgr_val)], temp_val, melt_total, color='blue')
#ax.scatter(temp_val, sgr_val, melt_total, color='blue')
ax.scatter(t, s, m, color='blue')
X, Y = np.meshgrid(t_range, s_range)
Z = func((X, Y), *popt)
ax.plot_surface(X, Y, Z, color='red', alpha=0.5)
ax.set_xlabel('temperature')
ax.set_ylabel('sgr')
ax.set_zlabel('melt')
plt.show()
'''

plt.figure(figsize=(4, 4))
clr = 'kbrcgy'
smb = 'o^spdh'
for t in range(len(temp)):
    plt.plot(sgr_val, melt_total[t,:], f'{clr[t]}:', linewidth=1, marker='o', fillstyle='none', markersize=5)
    plt.plot(sgr_val, meltTot_total[t,:], f'{clr[t]}--', linewidth=1, marker='x', markersize=5)
    z = k_fit*(temp_val[t] - To_fit) * s_range**(2/3)
    plt.plot(s_range, z, f'{clr[t]}-', linewidth=1, label = f'$T_b$ = {temp_val[t]}$^\circ$C')

plt.xlabel('$F_{s}$ (m$^3$/s)')
#plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.ylabel('Melt-flux anomaly (m$^3$/a)')
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
plt.tick_params(axis='both', labelsize=10)
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_L_vary_Fs_Tb_Fax_fit.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

plt.figure(figsize=(4, 4))
clr = 'kbrcgy'
smb = 'o^spdh'
for s in range(1, len(sgr)):
    plt.plot(temp_val, melt_total[:,s], f'{clr[s]}:', linewidth=1, marker='o', fillstyle='none', markersize=5)
    plt.plot(temp_val, meltTot_total[:,s], f'{clr[s]}--', linewidth=1 , marker='x', fillstyle='none', markersize=5)
    z = k_fit*(t_range - To_fit) * sgr_val[s]**(2/3)
    plt.plot(t_range, z, f'{clr[s]}-', linewidth=1, label = f'$F_s$ = {sgr_val_sd[s]}m$^3$/s')

plt.xlabel('$T_b$ ($^\circ$C)', fontsize=10)
#plt.ylabel('$\Delta \dot{m}$ (m/a)')
plt.ylabel('Melt-flux anomaly (m$^3$/a)', fontsize=10)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams['font.size'] = 10
plt.tick_params(axis='both', labelsize=10)
#plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_L_vary_Fs_Tb_Tax_fit.png', bbox_inches='tight', dpi=300)
else:
    plt.show()