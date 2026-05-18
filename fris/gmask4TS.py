#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import os

opt_save = 0
opt_noGL = 1

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

# MPAS Ocean data
fnum = 8

fris_loc = '/Users/ivankova/Desktop/Fris_hr'
f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'
mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'
out_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}melt_annual_2D_Y4.nc'

dsMesh = xarray.open_dataset(mesh_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell','maxLevelCell','layerThickness','cellsOnCell','nEdgesOnCell']]
dsMesh.load()

wct = dsMesh['layerThickness'].sum(dim='nVertLevels')
wct = np.squeeze(wct)

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
#iii = (FloatingMask == 1)

cellsOnCell = dsMesh.cellsOnCell.values
cellsOnCell = cellsOnCell

nEdgesOnCell = dsMesh.nEdgesOnCell.values
nEdgesOnCell = nEdgesOnCell

#Grounding line cells
nonzero_counts = np.count_nonzero(cellsOnCell, axis=1)
# iGL = np.where(nonzero_counts < nEdgesOnCell)[0]
iGL = (nonzero_counts < nEdgesOnCell)

print('Plotting')

fHeight = 5
fWidth = 6
plt.figure(figsize=(fWidth, fHeight))
#plt.plot(lon,lat,'b.')
#plt.plot(lon[iii],lat[iii],'r.')


#FRIS
iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
iam2 = (lon < 324.1) | (lat < -78.253)
iam = np.logical_and(iam1, iam2)
FiGL = iam & iGL
plt.plot(lon[FiGL], lat[FiGL], 'x', color='black')

if opt_noGL == 1:
    iam = iam & ~iGL

plt.plot(lon[iam], lat[iam], '.', color='lightgray')

print("Type iGL:", type(iGL))
bool_array = iam
print("Type FRIS:", type(iam))
print("Shape FRIS:", bool_array.shape)
count = np.sum(bool_array)
print("Number of cells where FRIS == 1:", count)


#FRIS Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (wct < 1500)
if opt_noGL == 1:
    iam = iam & ~iGL
plt.plot(lon[iam],lat[iam],'g.')
#Ronne Shelf
iam = (FloatingMask == 0) & (lat > -80) & (lat < -72.5) & (lon > 298) & (lon < 304.5) & (wct < 1200) & (wct > 500)
if opt_noGL == 1:
    iam = iam & ~iGL
plt.plot(lon[iam],lat[iam],'y.')
#Filchner Depression
iam = (FloatingMask == 0) & (lat > -80) & (lat < -74.9) & (lon > 316) & (lon < 332) & (wct < 1200) & (wct > 500)
if opt_noGL == 1:
    iam = iam & ~iGL
plt.plot(lon[iam],lat[iam],'m.')

#Berkner Bank
iam = (FloatingMask == 0) & (lat > -80) & (lat < -76) & (lon > 307) & (lon < 316) & (wct < 300)
if opt_noGL == 1:
    iam = iam & ~iGL
plt.plot(lon[iam],lat[iam],'r.')

#FRIS deep beneath
#iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.9) & (lon > 308) & (lon < 332) & (H < 2500) & (H > 1000)
iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
iam2 = (lon < 324.1) | (lat < -78.253)
iam3 = np.logical_and(iam1, iam2)
iam4 = (wct > 500)# & (H > 800)
iam = np.logical_and(iam3, iam4)
if opt_noGL == 1:
    iam = iam & ~iGL
plt.plot(lon[iam],lat[iam],'b.')

#DeepBack
iam = (FloatingMask == 1) & (lat > -84) & (lat < -80.5) & (lon > 300) & (lon < 316) & (wct > 500)# & (H > 800)
plt.plot(lon[iam],lat[iam],'y.')
if opt_noGL == 1:
    iam = iam & ~iGL

#DeepRonne
iam = (FloatingMask == 1) & (lat > -79) & (lat < -75.5) & (lon > 282) & (lon < 294) & (wct > 500)# & (H > 800)
plt.plot(lon[iam],lat[iam],'r.')
if opt_noGL == 1:
    iam = iam & ~iGL

#DeepFilchner
iam = (FloatingMask == 1) & (lat > -81.5) & (lat < -78) & (lon > 317) & (lon < 324) & (wct > 500)# & (H > 800)
plt.plot(lon[iam],lat[iam],'k.')
if opt_noGL == 1:
    iam = iam & ~iGL


'''
#FRIS Evans
iam = (FloatingMask == 1) & (lat > -77.2) & (lat < -75.75) & (lon > 281) & (lon < 285)
plt.plot(lon[iam],lat[iam],'y.')

#FRIS Rutford
iam = (FloatingMask == 1) & (lat > -79.6) & (lat < -78.2) & (lon > 275.5) & (lon < 281)
plt.plot(lon[iam],lat[iam],'b.')

#FRIS Institute
iam = (FloatingMask == 1) & (lat > -81) & (lat < -80.5) & (lon > 284) & (lon < 289)
plt.plot(lon[iam],lat[iam],'c.')

#FRIS Foundation
iam= (FloatingMask == 1) & (lat > -84) & (lat < -83.04) & (lon > 298) & (lon < 305)
plt.plot(lon[iam],lat[iam],'g.')

#FRIS Support Force
iam= (FloatingMask == 1) & (lat > -83) & (lat < -81.8) & (lon > 313.5) & (lon < 317)
plt.plot(lon[iam],lat[iam],'m.')

#FRIS Recovery
iam= (FloatingMask == 1) & (lat > -81) & (lat < -80.62) & (lon > 321) & (lon < 325)
plt.plot(lon[iam], lat[iam], '.', color='lightcoral')

#FRIS Slessor
iam= (FloatingMask == 1) & (lat > -80.5) & (lat < -79.9) & (lon > 327) & (lon < 332)
plt.plot(lon[iam],lat[iam],'r.')

#FRIS Bailey
iam= (FloatingMask == 1) & (lat > -79.7) & (lat < -79) & (lon > 326) & (lon < 332)
plt.plot(lon[iam],lat[iam],'k.')
'''

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tick_params(axis='both', labelsize=10)


if opt_save == 1:
    plt.savefig(f'{fris_loc}/Fris_plots/melt/Melt_hist_map_F{fnum}.png',
                bbox_inches='tight', dpi=300)
else:
    plt.show()



'''
#FRIS
iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
iam2 = (lon < 324.1) | (lat < -78.253)
iam = np.logical_and(iam1, iam2)

plt.plot(lon[iam],lat[iam],'k.')

#FRIS Evans
iam = (FloatingMask == 1) & (lat > -77.2) & (lat < -75.75) & (lon > 281) & (lon < 285)
plt.plot(lon[iam],lat[iam],'y.')

#FRIS Institute
iam = (FloatingMask == 1) & (lat > -81) & (lat < -79.5) & (lon > 283.5) & (lon < 289)
plt.plot(lon[iam],lat[iam],'c.')

#FRIS Rutford
iam = (FloatingMask == 1) & (lat > -79.6) & (lat < -78.2) & (lon > 275.5) & (lon < 281)
plt.plot(lon[iam],lat[iam],'b.')

#FRIS Foundation
iam= (FloatingMask == 1) & (lat > -84) & (lat < -82.3) & (lon > 298) & (lon < 305)
plt.plot(lon[iam],lat[iam],'g.')

#FRIS Support Force
iam= (FloatingMask == 1) & (lat > -83) & (lat < -81.5) & (lon > 313.5) & (lon < 317)
plt.plot(lon[iam],lat[iam],'m.')

#FRIS Recovery
iam= (FloatingMask == 1) & (lat > -81) & (lat < -80.5) & (lon > 321) & (lon < 325)
plt.plot(lon[iam],lat[iam],'g.')

#FRIS Slessor
iam= (FloatingMask == 1) & (lat > -80.3) & (lat < -79) & (lon > 324) & (lon < 332)
plt.plot(lon[iam],lat[iam],'r.')
'''


