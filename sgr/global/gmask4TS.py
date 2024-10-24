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

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'

#o_file = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_61-70_ts_1-70/mpaso_ANN_006101_007012_climo.nc'

dsMesh = xarray.open_dataset(p_file)
#dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','temperature','salinity','nVertLevels']]
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell','restingThickness']]
dsMesh.load()

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
iii = (FloatingMask == 1)

H = np.squeeze(dsMesh.restingThickness.data)
H = np.nansum(H, axis=1)
H = np.squeeze(H)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(lon,lat,'b.')
plt.plot(lon[iii],lat[iii],'r.')

#AMERY
iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
plt.plot(lon[iam],lat[iam],'g.')

#AMERY Shelf
iam = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 70) & (lon < 80) & (H < 1500)
plt.plot(lon[iam],lat[iam],'y.')

#PIG
iam = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880+360) & (lon < -98.4832+360)
plt.plot(lon[iam],lat[iam],'m.')
#print(len(lat[iam]))
#Thwaites
iam = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729+360) & (lon < -104.0573+360)
plt.plot(lon[iam],lat[iam],'c.')
#Amundsen shelf (Thwaites)
iam = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 256) & (lon < 260)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'y.')
#Amundsen shelf (PIG)
iam = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 252) & (lon < 256)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'k.')

#Getz
iam = (FloatingMask == 1) & (lat > -75) & (lat < -73.5) & (lon > 225) & (lon < 245.5)
plt.plot(lon[iam],lat[iam],'y.')
#Amundsen shelf (Getz upstream)
iam = (FloatingMask == 0) & (lat > -75) & (lat < -71) & (lon > 235) & (lon < 246)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'g.')


#Ross
iam = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212)
plt.plot(lon[iam],lat[iam],'m.')

#FRIS
iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
plt.plot(lon[iam],lat[iam],'k.')
#Ronne Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 315) & (H < 1500)
plt.plot(lon[iam],lat[iam],'y.')
#Filchner Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 316) & (lon < 332) & (H < 1500)
plt.plot(lon[iam],lat[iam],'m.')


#Totten
iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)
plt.plot(lon[iam],lat[iam],'k.')


plt.show()


