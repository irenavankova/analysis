#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

opt_save = 1
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

opt_23 = 1

t1 = 71
t1 = 61
t2 = t1+9
ts2 = t2

if opt_23 == 1:
    t1 = 23
    t2 = 26
    ts2 = 32


#trun = 'control'
#ttle = f'{trun} {t1}-{t2}'
#reg = 'Ross'
#reg = 'Amery'
#reg = 'PIG'
#reg = 'Thwaites'
reg = 'Totten'
ttle = f'{reg}_decade_{t1}_{t2}'
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
o_file_c = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{t1}-{t2}_ts_1-{ts2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'
o_file_m = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{t1}-{t2}_ts_1-{ts2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'

dsMesh = xarray.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell']]
dsMesh.load()
areaCell = np.squeeze(dsMesh.areaCell.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)
maxLevelCell = np.squeeze(dsMesh.maxLevelCell.data)-1

#Control
dsOutc = xarray.open_dataset(o_file_c)
dsOutc = dsOutc[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_landIceFreshwaterFluxTotal']]
dsOutc.load()
PTc = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_temperature.data)
PSc = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_salinity.data)
MRc = np.squeeze(dsOutc.timeMonthly_avg_landIceFreshwaterFluxTotal.data)

#MALI
dsOutm = xarray.open_dataset(o_file_m)
dsOutm = dsOutm[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_landIceFreshwaterFluxTotal']]
dsOutm.load()
PTm = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_temperature.data)
PSm = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_salinity.data) / rho_fw * secPerYear
MRm = np.squeeze(dsOutc.timeMonthly_avg_landIceFreshwaterFluxTotal.data) / rho_fw * secPerYear

#Get MASK
lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

reg = ["PIG", "Thwaites" , "Totten" , "Amery" , "Ross", "FRIS"]

for s in range(len(reg)):
    if reg[s] == 'Amery':
        iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
    elif reg[s] == 'PIG':
        iam = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880+360) & (lon < -98.4832+360)
    elif reg[s] == 'Thwaites':
        iam = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729+360) & (lon < -104.0573+360)
    elif reg[s] == 'Ross':
        iam = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212)
    elif reg[s] == 'FRIS':
        iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
    elif reg[s] == 'Totten':
        iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)

    mlc = maxLevelCell[iam]
    PTc_floor = np.zeros((len(mlc)))
    PTciam = PTc[iam,:]
    for v in range(len(mlc)):
        PTc_floor[v] = PTciam[v, mlc[v]]

    MFc = np.sum(MRc[iam] * areaCell[iam])
    MFm = np.sum(MRm[iam] * areaCell[iam])

    TFc = np.sum(PTc_floor * areaCell[iam])/np.sum(areaCell[iam])
    TFcmax = np.nanmax(PTc_floor)

    print(MFm-MFc)
    print(TFc)
    print(TFcmax)





