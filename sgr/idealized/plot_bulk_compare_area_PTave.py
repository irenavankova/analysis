#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

opt_save = 1
opt_Line = 1
plot_opt_ave = 0
hplot = 0

opt_thr = 3
t = 0
temp = ["rd", "rdn"]

d_scale = 30000

dir_fig_save = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/bulk/area/D{d_scale/1000}km_{temp[t]}_THR{opt_thr}'
if not os.path.exists(dir_fig_save):
    os.mkdir(dir_fig_save)

#h = 0
hloc = ["142", "122", "132"]
hloc_val = ["PE", "PC", "PW"]
#hloc = ["132"]
#hloc_val = ["PW"]
#hloc = ["142"]
#hloc_val = ["PE"]

ymaxA = np.array([0.2, 0.4, 0.8, 0.7])

if opt_Line == 1:
    hloc = ["112"]
    hloc_val = ["L"]
    plot_opt_ave = 0
    hplot = 0
    ymaxA[opt_thr] = 0.5
    ymaxA = np.array([0.5, 0.7, 0.8, 1])

do_lat = 0
v = 2

#------------------------------------------------------------------------------------------------
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

sgr = ["N", "R", "A" , "D" , "C" , "E" , "B"]
sgr_factor = 7.2169
sgr_val = np.array([0, 1, 10, 25 , 50 , 75 , 100])*sgr_factor

lat = ["0", "45", "65", "75", "85", "90"]
if do_lat == 1:
    sgr = ["N", "A", "D", "C", "B"]
    sgr_factor = 7.2169
    sgr_val = np.array([0, 10, 25, 50, 100]) * sgr_factor
else:
    if t == 0:
        v = 3
    else:
        v = 0

#--Get melt rate anomaly
melt_total = np.zeros((len(hloc),len(sgr)))
area_vel_anom = np.zeros((len(hloc),len(sgr)))
for h in range(len(hloc)):

    fdirf = f'{p_base}/rd/rd_{hloc[h]}R'
    df = xarray.open_dataset(f'{fdirf}/sgr_data.nc')
    df.load()
    subglacialRunoffFlux = np.squeeze(df.subglacialRunoffFlux.data)
    sgr_iii = subglacialRunoffFlux > 0
    sgr_pt = np.where(sgr_iii)[0]
    print(sgr_pt)

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

            if hloc[h] == '112':
                for j in range(len(sgr_pt)):
                    dist2point_now = np.absolute(
                        np.sqrt((x - x[sgr_pt[j]]) ** 2 + (y - y[sgr_pt[j]]) ** 2 + (z - z[sgr_pt[j]]) ** 2))
                    if j == 0:
                        dist2point = dist2point_now
                    else:
                        dist2point = np.minimum(dist2point, dist2point_now)

            else:
                dist2point = np.absolute(np.sqrt((x - x[sgr_pt]) ** 2 + (y - y[sgr_pt]) ** 2 + (z - z[sgr_pt]) ** 2))

            #dist2point = np.absolute(np.sqrt((x - x[sgr_pt]) ** 2 + (y - y[sgr_pt]) ** 2 + (z - z[sgr_pt]) ** 2))
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

        melt_total[h, s] = np.nansum(melt[iii] * area) / np.sum(area)
    melt_total[h, :] = melt_total[h, :] - melt_total[h, 0]

    if hloc_val[h] == "PC":
        area_PC = np.copy(area)
        melt_scat_PC = np.copy(melt_scat)
        vel_scat_PC = np.copy(vel_scat)
        iii_PC = np.copy(iii)
    elif hloc_val[h] == "PE":
        area_PE = np.copy(area)
        melt_scat_PE = np.copy(melt_scat)
        vel_scat_PE = np.copy(vel_scat)
        iii_PE = np.copy(iii)
    elif hloc_val[h] == "PW":
        area_PW = np.copy(area)
        melt_scat_PW = np.copy(melt_scat)
        vel_scat_PW = np.copy(vel_scat)
        iii_PW = np.copy(iii)
    elif hloc_val[h] == "L":
        area_L = np.copy(area)
        melt_scat_L = np.copy(melt_scat)
        vel_scat_L = np.copy(vel_scat)
        iii_L = np.copy(iii)

#--Define curve fitting functions

def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

def f_n_lin(x, a):
    return a * x

def f_n_lin_nori(x, a, b):
    return a * x + b

def f_n_pow_1_3_nori(x, a, b):
    return a * np.power(x, 1/3) + b

def f_n_pow_2_3_nori(x, a, b):
    return a * np.power(x, 2/3) + b

xfit = np.linspace(0,sgr_val[-1],100)

# PLOT and CALCULATE AREA WHERE MELT INCREASED BY MORE THAN A THRESHOLD
melt_anom_plume_mat = np.zeros((len(hloc), len(sgr) - 1))
melt_anom_out_mat = np.zeros((len(hloc), len(sgr) - 1))
area_anom_plume_mat = np.zeros((len(hloc), len(sgr) - 1))
area_anom_out_mat = np.zeros((len(hloc), len(sgr) - 1))
meltarea_anom_plume_mat = np.zeros((len(hloc), len(sgr) - 1))
meltarea_anom_out_mat = np.zeros((len(hloc), len(sgr) - 1))
area_tot_all = np.zeros(len(hloc))

for h in range(len(hloc)):
    if hloc_val[h] == "PC":
        area = np.copy(area_PC)
        melt_scat = np.copy(melt_scat_PC)
    elif hloc_val[h] == "PE":
        area = np.copy(area_PE)
        melt_scat = np.copy(melt_scat_PE)
    elif hloc_val[h] == "PW":
        area = np.copy(area_PW)
        melt_scat = np.copy(melt_scat_PW)
    elif hloc_val[h] == "L":
        area = np.copy(area_L)
        melt_scat = np.copy(melt_scat_L)

    area_anom = np.copy(area)
    area_tot_all[h] = np.sum(area)
    #plt.figure(figsize=(4, 4))
    for s in range(len(sgr)-1):
        melt_anom = np.copy(melt_scat[s, :])
        if opt_thr == 0:
            melt_thresh = np.max(melt_anom) / 10
        elif opt_thr == 1:
            melt_thresh = np.max(melt_anom) / 20
        elif opt_thr == 2:
            melt_thresh = 1
        elif opt_thr == 3:
            melt_thresh = 2
        #melt_thresh = np.max(melt_anom)
        #b_anom = np.absolute(melt_anom) >= melt_thresh
        b_anom_plume = melt_anom >= melt_thresh
        b_anom_out = np.copy(b_anom_plume)
        b_anom_out = np.invert(b_anom_out)
        #!!!!!
        #b_anom = melt_anom <= melt_thresh
        #plt.scatter(x_scat[b_anom_plume], y_scat[b_anom_plume], s=melt_anom[b_anom_plume] / np.max(np.max(melt_scat))*100, marker='o')
        #plt.scatter(x_scat[b_anom_out], y_scat[b_anom_out], s=melt_anom[b_anom_out] / np.max(np.max(melt_scat))*100, marker='o')

        melt_anom_plume_mat[h,s] = np.nansum(melt_anom[b_anom_plume] * area_anom[b_anom_plume]) / np.sum(area_anom[b_anom_plume])
        melt_anom_out_mat[h,s] = np.nansum(melt_anom[b_anom_out] * area_anom[b_anom_out]) / np.sum(area_anom[b_anom_out])
        area_anom_plume_mat[h,s] = np.sum(area_anom[b_anom_plume])
        area_anom_out_mat[h,s] = np.sum(area_anom[b_anom_out])
        meltarea_anom_plume_mat[h,s] = np.nansum(melt_anom[b_anom_plume] * area_anom[b_anom_plume])
        meltarea_anom_out_mat[h,s] = np.nansum(melt_anom[b_anom_out] * area_anom[b_anom_out])

        #plt.title(f'F = {sgr_val[s]}', fontsize = 8)
        #plt.show()
        #plt.show(block=False)

if plot_opt_ave == 0:
    h = hplot
    melt_anom_plume = np.copy(melt_anom_plume_mat[h,:])
    melt_anom_out = np.copy(melt_anom_out_mat[h,:])
    area_anom_plume = np.copy(area_anom_plume_mat[h,:])
    area_anom_out = np.copy(area_anom_out_mat[h,:])
    meltarea_anom_plume = np.copy(meltarea_anom_plume_mat[h,:])
    meltarea_anom_out = np.copy(meltarea_anom_out_mat[h,:])
    melt_tot = melt_total[h, :]
    area_tot = area_tot_all[h]
    ttl = hloc_val[h]
    ttl_save = hloc_val[h]
else:
    melt_anom_plume = np.mean(melt_anom_plume_mat, axis=0)
    melt_anom_out = np.mean(melt_anom_out_mat, axis=0)
    area_anom_plume = np.mean(area_anom_plume_mat, axis=0)
    area_anom_out = np.mean(area_anom_out_mat, axis=0)
    meltarea_anom_plume = np.mean(meltarea_anom_plume_mat, axis=0)
    meltarea_anom_out = np.mean(meltarea_anom_out_mat, axis=0)
    melt_tot = np.mean(melt_total, axis=0)
    area_tot = np.mean(area_tot_all)
    ttl = 'Channelized average'
    ttl_save = 'PTave'

#--PLOT Total melt rate flux and partition to "plume" and "ambient" (or high and low)
plt.figure(figsize=(4, 4))
clr = 'kbrkbc'
smb = 'o^o^^^'

ynow = melt_tot*area_tot
plt.plot(sgr_val, ynow, 'ko', linewidth=1, fillstyle='full', markersize=4, label = 'All')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow, p0=10**8)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'k--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow, p0=10**8)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'k-', linewidth=1)
print(popt)

ynow = np.insert(meltarea_anom_plume, 0, 0, axis=0)
plt.plot(sgr_val, ynow, 'ro', linewidth=1, fillstyle='full', markersize=4, label = 'High')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow, p0=10**8)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)

ynow = np.insert(meltarea_anom_out, 0, 0, axis=0)
plt.plot(sgr_val, ynow, 'bo', linewidth=1, fillstyle='full', markersize=4, label = 'Low')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'b-', linewidth=1)
popt, pcov = curve_fit(f_n_lin, sgr_val, ynow)
plt.plot(xfit, f_n_lin(xfit, *popt), 'b:', linewidth=1)

#plt.plot(tfit, qfit, 'k--', linewidth=1, label='fit')
plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('Area-integrated melt rate anomaly (m$^3$/a)')
plt.title(f'{ttl}, Lat = {lat[v]}S', fontsize = 8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.ylim([-0.1, 3.1])

if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_area_meltflux_{ttl_save}.png', bbox_inches='tight', dpi=300)
else:
    plt.show(block=False)

#--PLOT mean melt rate over an area that pertains to "plume" vs "ambient" (or high and low)

plt.figure(figsize=(4, 4))

ynow = melt_tot
plt.plot(sgr_val, ynow, 'ko', linewidth=1, fillstyle='full', markersize=4, label = 'All')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow, p0=10**8)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'k--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow, p0=10**8)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'k-', linewidth=1)

ynow = np.insert(melt_anom_plume, 0, 0, axis=0)
plt.plot(sgr_val, ynow, 'ro', linewidth=1, fillstyle='full', markersize=4, label = 'High')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'r-', linewidth=1)

ynow = np.insert(melt_anom_out, 0, 0, axis=0)
plt.plot(sgr_val, ynow, 'bo', linewidth=1, fillstyle='full', markersize=4, label = 'Low')
popt, pcov = curve_fit(f_n_pow_1_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, sgr_val, ynow)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'b-', linewidth=1)
popt, pcov = curve_fit(f_n_lin, sgr_val, ynow)
plt.plot(xfit, f_n_lin(xfit, *popt), 'b:', linewidth=1)

plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('Area-averaged melt rate anomaly (m/a)')
plt.title(f'{ttl}, Lat = {lat[v]}S', fontsize = 8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})

if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_area_meltave_{ttl_save}.png', bbox_inches='tight', dpi=300)
else:
    plt.show(block=False)


#--PLOT area fraction that pertains to "plume" vs "ambient" (or high and low)

plt.figure(figsize=(4, 4))

xnow = sgr_val[1:len(sgr_val)]

'''
ynow = area_anom_plume/np.sum(area) + area_anom_out/np.sum(area)
plt.plot(xnow, ynow, 'ko', linewidth=1, fillstyle='none', markersize=4)
'''

ynow = area_anom_plume/area_tot
plt.plot(xnow, ynow, 'ro', linewidth=1, fillstyle='full', markersize=4, label = 'High')
popt, pcov = curve_fit(f_n_pow_1_3_nori, xnow, ynow)
plt.plot(xfit, f_n_pow_1_3_nori(xfit, *popt), 'r--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3_nori, xnow, ynow)
plt.plot(xfit, f_n_pow_2_3_nori(xfit, *popt), 'r-', linewidth=1)
popt, pcov = curve_fit(f_n_lin_nori, xnow, ynow)
plt.plot(xfit, f_n_lin_nori(xfit, *popt), 'r:', linewidth=1)

'''
popt, pcov = curve_fit(f_n_pow_1_3, xnow, ynow)
plt.plot(xfit, f_n_pow_1_3(xfit, *popt), 'c--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3, xnow, ynow)
plt.plot(xfit, f_n_pow_2_3(xfit, *popt), 'c-', linewidth=1)
popt, pcov = curve_fit(f_n_lin, xnow, ynow)
plt.plot(xfit, f_n_lin(xfit, *popt), 'c:', linewidth=1)
'''

'''
ynow = area_anom_out/np.sum(area)
plt.plot(xnow, ynow, 'bo', linewidth=1, fillstyle='none', markersize=4)
popt, pcov = curve_fit(f_n_pow_1_3_nori, xnow, ynow)
plt.plot(xfit, f_n_pow_1_3_nori(xfit, *popt), 'b--', linewidth=1)
popt, pcov = curve_fit(f_n_pow_2_3_nori, xnow, ynow)
plt.plot(xfit, f_n_pow_2_3_nori(xfit, *popt), 'b-', linewidth=1)
popt, pcov = curve_fit(f_n_lin_nori, xnow, ynow)
plt.plot(xfit, f_n_lin_nori(xfit, *popt), 'b:', linewidth=1)
'''

plt.xlabel('$F_{s}$ (m$^3$/s)')
plt.ylabel('Area fraction')
plt.title(f'{ttl}, Lat = ${lat[v]}^\circ$S', fontsize = 8)
plt.legend(loc=2, prop={'size': 8})
plt.grid()
plt.rcParams.update({'font.size': 8})
plt.ylim([0, ymaxA[opt_thr]])

if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/plot_bulk_compare_area_areafrac_{ttl_save}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()

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

# MAKE SCATTER PLOT AND FIT A LINE

