#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 0

h = 0
hloc = ["122"]
hloc_val = ["PC"]
#hloc = ["132"]
#hloc_val = ["PW"]

d_scale = 30000

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

t = 1
temp = ["rd", "rdn"]
sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 25 , 50 , 75 , 100])*sgr_factor

do_lat = 0
v = 1
lat = ["0", "45", "65", "75", "85", "90"]
if do_lat == 1:
    sgr = ["N", "A", "D", "C", "B"]
    sgr_factor = 7.2169
    sgr_val = np.array([0, 10, 25, 50, 100]) * sgr_factor


fdirf = f'{p_base}/rd/rd_{hloc[h]}R'
df = xarray.open_dataset(f'{fdirf}/sgr_data.nc')
df.load()
subglacialRunoffFlux = np.squeeze(df.subglacialRunoffFlux.data)
sgr_iii = subglacialRunoffFlux > 0
sgr_pt = np.where(sgr_iii)[0]
print(sgr_pt)

melt_total = np.zeros((len(hloc),len(sgr)))
area_vel_anom = np.zeros((len(hloc),len(sgr)))
for h in range(len(hloc)):
    for s in range(len(sgr)):
        fdir = f'{p_base}/{temp[t]}/{temp[t]}_{hloc[h]}{sgr[s]}'
        if do_lat == 1:
            fdir = f'{p_base}/{temp[t]}L{lat[v]}/{temp[t]}/{temp[t]}_{hloc[h]}{sgr[s]}_L{lat[v]}'

        ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
        ds.load()
        melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
        #melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFlux.data)
        melt = melt / rho_fw * secPerYear
        vel = np.squeeze(ds.timeMonthly_avg_landIceFrictionVelocity.data)*100 #convert to cm/s

        if s == 0:
            dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
            dsMesh.load()
            areaCell = np.squeeze(dsMesh.areaCell.data)
            FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
            x = np.squeeze(dsMesh.xCell.data)
            y = np.squeeze(dsMesh.yCell.data)
            z = np.squeeze(dsMesh.zCell.data)
            dist2point = np.absolute(np.sqrt((x - x[sgr_pt]) ** 2 + (y - y[sgr_pt]) ** 2 + (z - z[sgr_pt]) ** 2))
            iii = (FloatingMask == 1) & (dist2point < d_scale)
            x_scat = x[iii]
            y_scat = y[iii]
            area = np.copy(areaCell[iii])
            melt_scat = np.zeros((len(sgr)-1, len(x[iii])))
            melt_ref = np.copy(melt[iii])
            vel_scat = np.zeros((len(sgr) - 1, len(x[iii])))
            vel_ref = np.copy(vel[iii])
            print('X runoff')
            print(x[sgr_pt])
            print('Y runoff')
            print(y[sgr_pt])
        else:
            melt_scat[s - 1,:] = np.copy(melt[iii]-melt_ref)
            vel_scat[s - 1, :] = np.copy(vel[iii] - vel_ref)

        #vel = vel * is_mask * 100 * i_ave_mask
        #melt = melt * is_mask * i_ave_mask
        #area = areaCell * is_mask * i_ave_mask
        melt = melt[iii]
        #print(np.shape(melt))
        #print(np.shape(area))
        melt_total[h, s] = np.nansum(melt * area) / np.sum(area)
    melt_total[h, :] = melt_total[h, :] - melt_total[h, 0]

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit = np.linspace(0,sgr_val[-1],100)

# FOR A GIVEN SGR, PLOT AREA OF MELT INCREASE above a threshold
s = len(sgr)-3
s = 0
clr = 'brkycmg'
smb = '^ox>sph'
melt_thresh_vec = np.linspace(0, 1, 11)
#melt_thresh_vec = np.array([0, 0.1, 0.9])
area_total = np.zeros(len(melt_thresh_vec)) + np.sum(area)
area_anom = np.zeros(len(area))
plt.figure(figsize=(4, 4))
for s in range(len(sgr)-1):
#for s in range(2):
    print('----------------------------------------------------------')
    print(sgr[s+1])
    area_thresh_vec = np.zeros(len(melt_thresh_vec))
    melt_anom = np.copy(melt_scat[s, :])
    print('MELT anom')
    # print(len(melt_anom))
    print(melt_anom)
    for m in range(len(melt_thresh_vec)):
        #print('MELT THRESH')
        #print(melt_thresh_vec[m]*np.max(melt_anom))
        b_anom = np.absolute(melt_anom) >= melt_thresh_vec[m]*np.max(melt_anom)
        #ind_x = np.where(np.absolute(melt_anom) < melt_thresh_vec[m]*np.max(melt_anom))
        #print('b_anom sum TRUE')
        #print(b_anom)
        #print('MELT anom')
        #print(len(melt_anom))
        #print(melt_anom)
        #print('ind_x LENGTH')
        #print(len(ind_x))
        area_anom = np.copy(area)
        area_thresh_vec[m] = np.sum(area_anom[b_anom])
        #print('area_anom[b_anom]')
        #print(area_anom[b_anom])
        #print('area[b_anom]')
        #print(area[b_anom])
        #print('AREA SUM ind_x')
        #print(np.sum(area_anom[b_anom]))
        #The following will introduce zeroes that can't get rid off!
        #area_anom[ind_x] = 0
        #print(np.sum(area_anom))
        #del b_anom, ind_x, area_anom
    plt.plot(melt_thresh_vec, area_thresh_vec, f'{clr[s]}:', linewidth=1, marker=f'{smb[s]}', fillstyle='none', markersize=4)

plt.plot(melt_thresh_vec, area_total, 'k--', linewidth=1, fillstyle='none', markersize=4)
plt.plot(melt_thresh_vec, area_total*0, 'k--', linewidth=1, fillstyle='none', markersize=4)
plt.xlabel('melt_thresh')
plt.ylabel('area_thresh')
plt.show(block=False)
#plt.show()


# PLOT AREA WHERE MELT INCREASED BY MORE THAN A THRESHOLD
melt_anom_plume = np.zeros(len(sgr)-1)
melt_anom_out = np.zeros(len(sgr)-1)
area_anom_plume = np.zeros(len(sgr)-1)
area_anom_out = np.zeros(len(sgr)-1)
meltarea_anom_plume = np.zeros(len(sgr)-1)
meltarea_anom_out = np.zeros(len(sgr)-1)
area_anom = np.copy(area)
plt.figure(figsize=(4, 4))
for s in range(len(sgr)-1):
    melt_anom = np.copy(melt_scat[s, :])
    melt_thresh = np.max(melt_anom)/20
    #melt_thresh = np.max(melt_anom)
    #melt_thresh = 1.0
    #b_anom = np.absolute(melt_anom) >= melt_thresh
    b_anom_plume = melt_anom >= melt_thresh
    b_anom_out = np.copy(b_anom_plume)
    b_anom_out = np.invert(b_anom_out)
    #!!!!!
    #b_anom = melt_anom <= melt_thresh
    plt.scatter(x_scat[b_anom_plume], y_scat[b_anom_plume], s=melt_anom[b_anom_plume] / np.max(np.max(melt_scat))*100, marker='o')
    plt.scatter(x_scat[b_anom_out], y_scat[b_anom_out], s=melt_anom[b_anom_out] / np.max(np.max(melt_scat))*100, marker='o')

    melt_anom_plume[s] = np.nansum(melt_anom[b_anom_plume] * area_anom[b_anom_plume]) / np.sum(area_anom[b_anom_plume])
    melt_anom_out[s] = np.nansum(melt_anom[b_anom_out] * area_anom[b_anom_out]) / np.sum(area_anom[b_anom_out])
    area_anom_plume[s] = np.sum(area_anom[b_anom_plume])
    area_anom_out[s] = np.sum(area_anom[b_anom_out])
    meltarea_anom_plume[s] = np.nansum(melt_anom[b_anom_plume] * area_anom[b_anom_plume])
    meltarea_anom_out[s] = np.nansum(melt_anom[b_anom_out] * area_anom[b_anom_out])

    plt.title(f'F = {sgr_val[s]}', fontsize = 8)
    #plt.show()
    plt.show(block=False)

plt.figure(figsize=(4, 4))
plt.plot(sgr_val[1:len(sgr_val)], melt_anom_plume, 'ro', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], melt_anom_plume)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], melt_anom_plume)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)
plt.plot(sgr_val[1:len(sgr_val)], melt_anom_out, 'bo', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], melt_anom_out)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], melt_anom_out)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'b-', linewidth=1)
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
#plt.show()
plt.show(block=False)

plt.figure(figsize=(4, 4))
plt.plot(sgr_val[1:len(sgr_val)], meltarea_anom_plume, 'ro', linewidth=1, fillstyle='none', markersize=4)
plt.plot(sgr_val[1:len(sgr_val)], melt_anom_plume*area_anom_plume, 'kx', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], meltarea_anom_plume)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], meltarea_anom_plume)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)
plt.plot(sgr_val[1:len(sgr_val)], meltarea_anom_out, 'bo', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], meltarea_anom_out)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], meltarea_anom_out)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'b-', linewidth=1)
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a) x Area')
plt.show(block=False)

plt.figure(figsize=(4, 4))
plt.plot(sgr_val[1:len(sgr_val)], area_anom_plume, 'ro', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], area_anom_plume)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], area_anom_plume)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)
plt.plot(sgr_val[1:len(sgr_val)], area_anom_out, 'bo', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val[1:len(sgr_val)], area_anom_out)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val[1:len(sgr_val)], area_anom_out)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'b-', linewidth=1)
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('Area (m$^2$)')
plt.show(block=False)

'''
def fit_vel(x, a):
    return a * np.power(x-459000,1) + 39301.27018922

xfit = np.linspace(458000,500000,100)

#plt.figure(figsize=(4, 4))
for s in range(len(sgr)-1):
    plt.figure(figsize=(4, 4))

    vel_anom = vel_scat[s, :]
    i_anom = vel_anom < 0.02
    vel_anom[i_anom] = 0.001
    w = np.power(1/np.absolute(vel_anom),2)

    #print(i_anom)

    #plt.scatter(x_scat, y_scat, s=melt_scat[s,:], marker='o')
    plt.scatter(x_scat, y_scat, s=vel_anom/np.max(np.max(vel_scat)), marker='o')
    #plt.scatter(x_scat[i_anom], y_scat[i_anom], s=vel_scat[s, i_anom]/np.max(np.max(vel_scat)), marker='o')
    popt, pcov = curve_fit(fit_vel, x_scat, y_scat, sigma=w)
    plt.plot(xfit, fit_vel(xfit, *popt), 'k-', label='fit: a=%5.3f' % tuple(popt))
    #popt, pcov = curve_fit(fit_vel, x_scat, y_scat, sigma=np.power(np.absolute(vel_scat[s, :]),2))
    #plt.plot(xfit, fit_vel(xfit, *popt), 'r-', label='fit: a=%5.3f' % tuple(popt))
    plt.ylim([np.min(np.min(y_scat)), np.max(np.max(y_scat))])
    plt.show()
#plt.show()
'''


# PLOT MEAN MELT RATE OVER A RESTRICTED AREA AROUND PLUME ORIGIN
plt.figure(figsize=(4, 4))
clr = 'brkbc'
smb = '^o^^^'
h = 0
plt.plot(sgr_val, melt_total[h,:], f'{clr[h]}:', linewidth=1, marker=f'{smb[h]}', fillstyle='none', markersize=4, label = f'{hloc_val[h]}')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, melt_total[h,:])
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, melt_total[h,:])
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('$\Delta \dot{m}$ (m/a)')
if temp[t] == "rd":
    plt.title('$T_b=1^\circ$C, $f=-1.409 \cdot 10^{-4}$ s$^{-1}$', fontsize = 8)
else:
    plt.title('$T_b=1^\circ$C, $f=0$ s$^{-1}$', fontsize=8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_horizontal_sgr_{temp[t]}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

# MAKE SCATTER PLOT AND FIT A LINE

