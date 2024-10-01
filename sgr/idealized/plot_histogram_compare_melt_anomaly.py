#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

opt_save = 0

h = 0
#hloc = ["122"]
#hloc_val = ["PC"]
#hloc = ["132"]
#hloc_val = ["PW"]
hloc = ["142"]
hloc_val = ["PE"]

t = 1
temp = ["rd", "rdn"]

do_lat = 0
v = 2

d_scale = 30000

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

fdirf = f'{p_base}/rd/rd_{hloc[h]}R'
df = xarray.open_dataset(f'{fdirf}/sgr_data.nc')
df.load()
subglacialRunoffFlux = np.squeeze(df.subglacialRunoffFlux.data)
sgr_iii = subglacialRunoffFlux > 0
sgr_pt = np.where(sgr_iii)[0]
print(sgr_pt)

#--Get melt rate anomaly
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

# PLOT and CALCULATE AREA WHERE MELT INCREASED BY MORE THAN A THRESHOLD
area_anom = np.copy(area)
plt.figure(figsize=(4, 4))
for s in range(len(sgr)-1):
    melt_anom = np.copy(melt_scat[s, :])

    plt.hist(melt_anom, density=True, bins=30)  # density=False would make counts

plt.ylabel('Probability')
plt.xlabel('melt');

#plt.title(f'F = {sgr_val[s]}', fontsize = 8)
plt.show()
#plt.show(block=False)

# MAKE SCATTER PLOT AND FIT A LINE

