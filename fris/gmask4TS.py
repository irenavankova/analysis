#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

opt_save = 0
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

# MPAS Ocean data
p_file = f'/Users/irenavankova/Desktop/Fris_ncfiles/F8mesh.nc'

dsMesh = xarray.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell','maxLevelCell']]
dsMesh.load()

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
iii = (FloatingMask == 1)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(lon,lat,'b.')
plt.plot(lon[iii],lat[iii],'r.')


#FRIS
iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)

iam2 = (lon < 324.1) | (lat < -78.253)
iam = np.logical_and(iam1, iam2)

plt.plot(lon[iam],lat[iam],'k.')


plt.show()


