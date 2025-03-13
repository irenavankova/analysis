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

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(lon,lat,'b.')
plt.plot(lon[iii],lat[iii],'r.')

#ANT-AMERY
iam = (FloatingMask == 1) & ((lat < -74.3608) | (lat > -67.7122) | (lon < 62.0419) | (lon > 78))
#plt.plot(lon[iam],lat[iam],'g.')
#plt.show()

#AMERY
iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
plt.plot(lon[iam],lat[iam],'g.')

#AMERY Shelf
iam = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 67.5) & (lon < 80) & (H < 1500)
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
#iam = (FloatingMask == 0) & (lat > -75) & (lat < -71) & (lon > 235) & (lon < 246)  & (H < 1500)
#Amundsen shelf (Getz all)
iam = (FloatingMask == 0) & (lat > -75.5) & (lat < -71) & (lon > 225) & (lon < 252)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'g.')

# Amundsen sea shelf all
iam = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 225) & (lon < 260)  & (H < 1500)
plt.plot(lon[iam],lat[iam],'.',color='gray')
# Amundsen sea ice shelves all
iam = (FloatingMask == 1) & (lat > -76) & (lat < -73.2) & (lon > 225) & (lon < 261)
plt.plot(lon[iam],lat[iam],'.',color='black')




#Ross
iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
iam2 = (lat < -77.8) | (lon < 200)
iam = np.logical_and(iam1, iam2)
plt.plot(lon[iam],lat[iam],'m.')
#Ross shelf
iam = (FloatingMask == 0) & (lat > -85.6) & (lat < -73) & (lon > 158.64) & (lon < 210) & (H < 1500)
plt.plot(lon[iam],lat[iam],'k.')


#FRIS
iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
plt.plot(lon[iam],lat[iam],'k.')
#Ronne Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 315) & (H < 1500)
plt.plot(lon[iam],lat[iam],'y.')
#Filchner Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 316) & (lon < 332) & (H < 1500)
plt.plot(lon[iam],lat[iam],'m.')
#FRIS Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (H < 1500)
plt.plot(lon[iam],lat[iam],'g.')

##Larsen
iam = (FloatingMask == 1) & (lat > -69.5) & (lat < -66.1) & (lon > 294) & (lon < 300)
plt.plot(lon[iam],lat[iam],'k.')

##Larsen shelf
iam = (FloatingMask == 0) & (lat > -69.5) & (lat < -66.1) & (lon > 298) & (lon < 310) & (H < 1500)
plt.plot(lon[iam],lat[iam],'m.')


##GeorgeVI
iam = (FloatingMask == 1) & ((lat > -74) & (lat < -70) & (lon > 287.3) & (lon < 293.4)) & ((lat < -72.5) | (lon > 290))
plt.plot(lon[iam],lat[iam],'k.')

##Stange
iam = (FloatingMask == 1) & (lat > -73.5) & (lat < -72.5) & (lon > 271) & (lon < 274)
plt.plot(lon[iam],lat[iam],'y.')

##Abbot
iam = (FloatingMask == 1) & ((lat > -73.34) & (lat < -2) & (lon > 256) & (lon < 271)) & ((lat < -72.28) | (lon < 260)) & ((lat > -73.2) | (lon > 260))
plt.plot(lon[iam],lat[iam],'c.')




#Totten
iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)
plt.plot(lon[iam],lat[iam],'k.')

#MoscowU
iam = (FloatingMask == 1) & ((lat > -68) & (lat < -66) & (lon > 119.5) & (lon < 122.3)) & ((lat < -67) | (lon > 120.5))
plt.plot(lon[iam],lat[iam],'k.')
#Totten + Moscow shelf
iam = (FloatingMask == 0) & (lat > -68) & (lat < -65) & (lon > 114) & (lon < 124) & (H < 1500)
plt.plot(lon[iam],lat[iam],'m.')
#Totten + Moscow
iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 122.3)
plt.plot(lon[iam],lat[iam],'g.')



#FImbul
iam = (FloatingMask == 1) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8))
plt.plot(lon[iam],lat[iam],'k.')
#FImbul shelf
iam = (FloatingMask == 0) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8)) & (H < 1500)
plt.plot(lon[iam],lat[iam],'m.')

#Nivlisen
iam = (FloatingMask == 1) & (lat > -71) & (lat < -69.8) & (lon > 9.5) & (lon < 12.9)
plt.plot(lon[iam],lat[iam],'c.')

#Roi Baudouin
iam = (FloatingMask == 1) & (lat > -71.5) & (lat < -69) & (lon >24) & (lon < 33)
plt.plot(lon[iam],lat[iam],'y.')

#Muninisen (Borchgrevink)
iam = (FloatingMask == 1) & (lat > -71.5) & (lat < -69) & (lon >19) & (lon < 22)
plt.plot(lon[iam],lat[iam],'m.')

#Ekstrom
iam = (FloatingMask == 1) & (lat > -72) & (lat < -70) & (lon >350) & (lon < 352.4)
plt.plot(lon[iam],lat[iam],'m.')

#Riiser-Larsen/Brunt
iam = (FloatingMask == 1) & (lat > -76) & (lat < -71.66) & (lon >332.5) & (lon < 350)
plt.plot(lon[iam],lat[iam],'c.')

#Shackleton
iam = (FloatingMask == 1) & (lat > -67) & (lat < -64) & (lon >92.5) & (lon < 105)
plt.plot(lon[iam],lat[iam],'c.')

#West
iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon >80) & (lon < 90)
plt.plot(lon[iam],lat[iam],'k.')

#Nickerson/Sulzberger/Swinbourne
iam = (FloatingMask == 1) & (lat > -78) & (lat < -75) & (lon > 206) & (lon < 220)
plt.plot(lon[iam],lat[iam],'c.')



plt.show()


