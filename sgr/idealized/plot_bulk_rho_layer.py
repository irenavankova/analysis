#!/usr/bin/env python3
import numpy as np
import xarray
import os
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'
sgr_base = '/Users/irenavankova/Work/data_sim/Compass_files/mesh_2k/sgr'

temp = 'rd'
sgr = ["R", "A", "A", "A", "A"]
#sgr = ["R", "B", "B", "B", "B"]
hloc = ["112", "112", "132", "122", "142"]

for c in range(len(sgr)):
    #Load files
    fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'

    dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
    dsMesh.load()
    minLevelCell = np.squeeze(dsMesh.minLevelCell.data)
    maxLevelCell = np.squeeze(dsMesh.maxLevelCell.data)
    nVert = np.squeeze(dsMesh.nVertLevels.data)
    areaCell = np.squeeze(dsMesh.areaCell.data)
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
    iii = FloatingMask == 1

    ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
    ds.load()
    rho = np.squeeze(ds.timeMonthly_avg_potentialDensity.data)
    if c == 0:
        rho_layer = np.zeros((len(sgr), len(nVert)))
        rho_layer_anom = np.zeros((len(sgr), len(nVert)))

    for v in range(len(nVert)):
        rho_layer[c, v] = np.nansum(FloatingMask[iii] * rho[iii,v] * areaCell[iii]) / np.sum(areaCell[iii])
    rho_layer_anom[c, :] = rho_layer[c, :] - rho_layer[0, :]

plt.figure(figsize=(4, 4))
for c in range(1, len(sgr)):
    plt.plot(nVert, rho_layer_anom[c, :], '-', label = f'{hloc[c]}{sgr[c]}')

plt.legend(loc=2, prop={'size': 6})
plt.show()
