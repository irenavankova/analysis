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
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
dsMesh = xarray.open_dataset(p_file)
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

# Ship CTD data: Shenjie Zhou
sz_file = f'/Users/irenavankova/Desktop/Shenjie/Shenjie_CT.nc'
dsSZ = xarray.open_dataset(sz_file)
dsSZ = dsSZ[['lat', 'lon','bathy']]
dsSZ.load()
latsz = np.squeeze(dsSZ.lat.data)
lonsz = np.squeeze(dsSZ.lon.data)
lonsz[lonsz < 0] = lonsz[lonsz < 0] + 360
Hsz = np.squeeze(dsSZ.bathy.data)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(lonsz,latsz,'b.')

#AMERY Shelf
#iam = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 67.5) & (lon < 80) & (H < 1500)
iam = (latsz > -70) & (latsz < -65) & (lonsz > 67.5) & (lonsz < 80) & (-Hsz < 1500)
print(iam.shape)
plt.plot(lonsz[iam],latsz[iam],'y.')

#Ross Shelf
iam = (latsz > -85.6) & (latsz < -73) & (lonsz > 158.64) & (lonsz < 210) & (-Hsz < 1500)
print(iam.shape)
plt.plot(lonsz[iam],latsz[iam],'k.')

#FRIS Shelf
iam = (latsz > -80) & (latsz < -72) & (lonsz > 298) & (lonsz < 332) & (-Hsz < 1500)
print(iam.shape)
plt.plot(lonsz[iam],latsz[iam],'g.')

#Amundsen Shelf
iam = (latsz > -76) & (latsz < -71) & (lonsz > 225) & (lonsz < 260)  & (-Hsz < 1500)
plt.plot(lonsz[iam],latsz[iam],'.',color='gray')



plt.show()


