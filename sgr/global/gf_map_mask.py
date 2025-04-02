#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.path as mpath

import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import cmocean
import gmask_is

opt_save = 1

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
dsMesh = xr.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell','restingThickness']]
dsMesh.load()

FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
lat = dsMesh['latCell'].values
lon = dsMesh['lonCell'].values
lat = lat*180/np.pi
lon = lon*180/np.pi
regMask = np.zeros([len(lat)])

y1 = '41'
y2 = '50'
fdir = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{y1}-{y2}_ts_1-{y2}/fw'
p_file = f'{fdir}/mpaso_ANN_{str(y1).zfill(4)}01_{str(y2).zfill(4)}12_climo.nc'

ds = xr.open_dataset(p_file)
layerThicknessC = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)
sshC = np.squeeze(ds['timeMonthly_avg_ssh'].values)
depthC = layerThicknessC
depthC[:,0] = depthC[:,0] - sshC
depthC = np.cumsum(depthC, axis=1)
depthC = np.nanmax(depthC, axis=1)

icav = (FloatingMask == 1)
ishelf = ((depthC < 1500) & (lat < -60.0)) | (FloatingMask == 1)
is60 = (lat < -60.0)

ctr = 0
#plt.plot(lon, lat, 'b.')

ctr = ctr + 1
regMask[is60] = ctr
#plt.plot(lon[is60], lat[is60], '.')

ctr = ctr + 1
regMask[ishelf] = ctr
#plt.plot(lon[ishelf], lat[ishelf], '.')

ctr = ctr + 1
regMask[icav] = ctr
#plt.plot(lon[icav], lat[icav], '.')

region_name = ["Amundsen","Amery","Ross","FRIS","Larsen","Fimbul"]
#region_name = ["FRIS"]

for j in range(len(region_name)):
    if region_name[j] == "Amery":
        ttl = "Amery"
        iceshelves = ["Amery_shelf", "Amery"]
        nshelf = 0
    elif region_name[j] == "Ross":
        ttl = "Ross"
        iceshelves = ["Ross_shelf", "Ross"]
        nshelf = 0
    elif region_name[j] == "FRIS":
        ttl = "Filchner-Ronne"
        iceshelves = ["Filchner-Ronne_shelf", "Filchner-Ronne"]
        nshelf = 0
    elif region_name[j] == "Amundsen":
        ttl = "Amundsen Sea shelves"
        iceshelves = ["Amundsen_shelf", "Amundsen"]
        nshelf = 0
    elif region_name[j] == "Larsen":
        ttl = "Larsen C"
        iceshelves = ["Larsen_C_shelf", "Larsen_C"]
        nshelf = 0
    elif region_name[j] == "Fimbul":
        ttl = "Fimbul"
        iceshelves = ["Fimbul_shelf", "Fimbul"]
        nshelf = 0

    iam, areaCell, isz = gmask_is.get_mask(iceshelves)

    for n in range(len(iceshelves)):
        ctr = ctr + 1
        iis = iam[n,:]
        regMask[iis] = ctr

        #plt.plot(lon[iis], lat[iis], '.')


#plt.show()

colors = ["lightgray", "slategray", "cyan","orange","yellow","brown","lightcoral","green","yellowgreen","royalblue","lightskyblue","darkorchid","plum","olive","darkkhaki"]

# Create the colormap
cmap_name = "segcor"
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors)


inum = regMask > 0

ig = plt.figure(num=None, figsize=(8, 6), edgecolor='k')

ax = plt.axes(projection=ccrs.SouthPolarStereo())

ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())

theta = np.linspace(0, 2*np.pi, 200)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

ax.set_boundary(circle, transform=ax.transAxes)

#ax.imshow(data.T, origin='lower', extent=[-180,180,-90,90], transform=ccrs.PlateCarree(),cmap='jet',vmin=0, vmax=1.0)
sc = ax.scatter(lon[inum], lat[inum], c=regMask[inum], cmap=custom_cmap, marker='.', s=1, edgecolor='none', transform=ccrs.PlateCarree())


#fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': ccrs.SouthPolarStereo()})
#ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())
#sc = ax.scatter(lon[inum], lat[inum], c=regMask[inum], cmap=custom_cmap, s=0.1, transform=ccrs.PlateCarree())

# Add a colorbar for the field A
#plt.colorbar(sc, ax=ax, label='Field A')
gridlines = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', linestyle='--')

# Label both latitude and longitude gridlines
#gridlines.xlabels_top = True  # Disable top longitude labels
#gridlines.xlabels_bottom = True  # Enable bottom longitude labels
#gridlines.ylabels_left = True  # Enable left latitude labels
#gridlines.ylabels_right = True  # Disable right latitude labels
#gridlines.xlabels_left = True  # Disable top longitude labels
#gridlines.xlabels_right = True  # Disable top longitude labels
gridlines.top_labels = True  # Add labels on the top
gridlines.bottom_labels = True # Add labels on the bottom
gridlines.left_labels = True # Add labels on the left
gridlines.right_labels = True # Add labels on the right
fsize = 8
gridlines.xlabel_style = {'size': fsize, 'color': 'gray', 'weight': 'normal'}
gridlines.ylabel_style = {'size':fsize, 'color': 'gray', 'weight': 'normal'}

# Show the plot

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/map/TS_regions.png', bbox_inches='tight', dpi=600)
else:
    plt.show()