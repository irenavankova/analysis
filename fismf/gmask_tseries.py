#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

sec = 'sea_ice'

opt_save = 0
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

# MPAS Ocean data
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
dsMesh = xarray.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell','restingThickness']]
dsMesh.load()

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
#plt.plot(lon,lat,'b.')

FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

if sec == 'ice_shelves':
    H = np.squeeze(dsMesh.restingThickness.data)
    H = np.nansum(H, axis=1)
    H = np.squeeze(H)
    iii = (FloatingMask == 1)
    plt.plot(lon[iii],lat[iii],'r.')


# Sea-ice regions

#Southern Ocean
iam = (lat < -50) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'b.')

#East Ant
iam = (lat < -60) & (lon > 90) & (lon < 150) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'c.')

#Amery
iam = (lat < -60) & (lon > 60) & (lon < 90) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'r.')

#DML
iam1 = (lat < -60) & (lon > 350) & (FloatingMask == 0)
iam2 = (lat < -60) & (lon > 0) & (lon < 60) & (FloatingMask == 0)
iam = np.logical_or(iam1, iam2)
plt.plot(lon[iam],lat[iam],'g.')

#Weddell
iam = (lat < -60) & (lon > 297) & (lon < 350) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'r.')

#AB
iam = (lat < -68) & (lon > 220) & (lon < 294) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'y.')

#Ross West
iam = (lat < -60) & (lon > 180) & (lon < 220) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'m.')

#Ross East
iam = (lat < -60) & (lon > 150) & (lon < 180) & (FloatingMask == 0)
plt.plot(lon[iam],lat[iam],'k.')

'''
# Ice shelf regions
#Peninnsula West -Bellingshausen
iam1 = (FloatingMask == 1) & (lat > -77.4) & (lon > 261) & (lon < 293.8)
iam2 = (FloatingMask == 1) & (lat > -77.3) & (lon > 255) & (lon < 262.8)
iam = np.logical_or(iam1, iam2)
plt.plot(lon[iam],lat[iam],'m.')

#Amundsen (up to Getz)
iam = (FloatingMask == 1) & (lat > -76) & (lat < -73.2) & (lon > 225) & (lon < 261)
plt.plot(lon[iam],lat[iam],'c.')

#Ross
iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
iam2 = (lat < -77.8) | (lon < 200)
iam = np.logical_and(iam1, iam2)
plt.plot(lon[iam],lat[iam],'m.')

#East Antarctica up to Amery
iam = (FloatingMask == 1) & (lat > -72) & (lon > 80) & (lon < 166)
plt.plot(lon[iam],lat[iam],'c.')

#Amery
iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
plt.plot(lon[iam],lat[iam],'g.')

#DML
iam = (FloatingMask == 1) & ((lon < 62.) | (lon > 332))
plt.plot(lon[iam],lat[iam],'y.')

#FRIS
iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
plt.plot(lon[iam],lat[iam],'k.')

#Peninnsula East -Larsens
iam = (FloatingMask == 1) &  (lat > -74.3) & (lon > 294) & (lon < 301)
plt.plot(lon[iam],lat[iam],'c.')
'''




'''
# Continental Shelves

# Amundsen sea shelf
iam = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 225) & (lon < 260)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'.',color='gray')

#Ross shelf
iam = (FloatingMask == 0) & (lat > -85.6) & (lat < -73) & (lon > 158.64) & (lon < 210) & (H < 1500)
plt.plot(lon[iam],lat[iam],'k.')

#AMERY Shelf
iam = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 67.5) & (lon < 80) & (H < 1500)
plt.plot(lon[iam],lat[iam],'y.')

#FRIS Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (H < 1500)
plt.plot(lon[iam],lat[iam],'g.')

'''



plt.show()


