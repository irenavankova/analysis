#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = 'rd'
sgr = ["A" , "D" , "C" , "E" , "B"]
hloc = '112'

for s in range(len(sgr)):
    fdir = f'{p_base}/{temp}/{temp}_{hloc}{sgr[s]}'

    ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
    #ds = ds['timeMonthly_avg_landIceFreshwaterFluxTotal']
    ds.load()
    melt = np.squeeze(ds.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
    melt = melt/rho_fw*secPerYear

    dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
    #dsMesh = dsMesh['nCells', 'landIceFloatingMask']
    dsMesh.load()
    #nCells = dsMesh.nCells.data
    areaCell = np.squeeze(dsMesh.areaCell.data)
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

    iii = FloatingMask == 1
    print(melt.shape)
    print(FloatingMask.shape)
    print(areaCell.shape)

    melt_total = np.sum(FloatingMask[iii]*melt[iii]*areaCell[iii])/np.sum(areaCell[iii])
    print(melt_total)
    #print(nCells)





