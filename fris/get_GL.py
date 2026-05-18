#!/usr/bin/env python3

import xarray as xr
import numpy as np
import gmask_reg
import os

import matplotlib.pyplot as plt

mesh_file = '/Users/ivankova/Desktop/Fris_hr/E3SM_init/ocean.ECwISC30to60E2r1.230220.nc'

dsMesh = xr.open_dataset(mesh_file)
#dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell','maxLevelCell','restingThickness','layerThickness']]
#dsMesh.load()

wct = dsMesh['layerThickness'].sum(dim='nVertLevels')
wct = np.squeeze(wct)

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

fHeight = 5
fWidth = 6
plt.figure(figsize=(fWidth, fHeight))

#FRIS
iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
iam2 = (lon < 324.1) | (lat < -78.253)
iam = np.logical_and(iam1, iam2)

plt.plot(lon[iam], lat[iam], '.', color='lightgray')

cellsOnCell = dsMesh.cellsOnCell.values
cellsOnCell = cellsOnCell

nEdgesOnCell = dsMesh.nEdgesOnCell.values
nEdgesOnCell = nEdgesOnCell

nonzero_counts = np.count_nonzero(cellsOnCell, axis=1)
#iGL = np.where(nonzero_counts < nEdgesOnCell)[0]
iGL = (nonzero_counts < nEdgesOnCell)

'''
print(cellsOnCell.shape)
for j, cell in enumerate(cellsOnCell):
    print(nEdgesOnCell[j])
    print(nonzero_counts[j])
    print(cell)
'''

#FRIS
iam = np.logical_and(iam, iGL)
#iam = np.logical_and(iam1, iam2, iGL)

plt.plot(lon[iam], lat[iam], '.', color='black')

plt.show()

