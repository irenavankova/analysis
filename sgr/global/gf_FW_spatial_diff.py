#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import cmocean
y1 = '50'
y2 = '50'
y1 = '110'
y2 = '110'

#opt_nlay40 = 40
#nlayleg = '_upper1000m'
opt_nlay40 = 0
nlayleg = '_all'

opt_mali = 1

if opt_mali == 8:
    mali_file = 'S12_mali_x8'
    fname = 'M8'
    cmax = 0.01
if opt_mali == 4:
    mali_file = 'S12_mali_x4'
    fname = 'M4'
    cmax = 0.08
if opt_mali == 1:
    mali_file = 'S12_mali'
    fname = 'M1'
    cmax = 0.01


# FW Ant projected file
fdir = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{y1}-{y2}_ts_1-{y2}/fw'
p_file = f'{fdir}/fw_9000.0x9000.0km_15.0km_Antarctic_stereo.nc'
ds = xr.open_dataset(p_file)
lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
srfwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values)
FWC = lifwc + srfwc
#FWC = lifwc
layerThicknessC = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)
sshC = np.squeeze(ds['timeMonthly_avg_ssh'].values)
depthC = layerThicknessC
depthC[:,:,0] = depthC[:,:,0] - sshC
depthC = np.cumsum(depthC, axis=2)
FWC[FWC < 0] = 0.0
FWC = np.nan_to_num(FWC, nan=0.0)

fdir = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{mali_file}/clim_{y1}-{y2}_ts_1-{y2}/fw'
p_file = f'{fdir}/fw_9000.0x9000.0km_15.0km_Antarctic_stereo.nc'
ds = xr.open_dataset(p_file)
lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
srfwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values)
FWM = lifwc + srfwc
#FWM = lifwc
layerThicknessM = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)
sshM = np.squeeze(ds['timeMonthly_avg_ssh'].values)
depthM = layerThicknessM
depthM[:,:,0] = depthM[:,:,0] - sshM
depthM = np.cumsum(depthM, axis=2)
FWM[FWC < 0] = 0.0
FWM = np.nan_to_num(FWM, nan=0.0)

mkm = 1000
x = ds['x'].values/mkm
y = ds['y'].values/mkm
lat = ds['lat'].values

xmax = 3500
xmin = -xmax
ymax = xmax
ymin = -ymax


if opt_nlay40 == 40:
    FWC = FWC[:, :, 0:39]
    FWM = FWM[:, :, 0:39]
    depthC = depthC[:, :, 0:39]
    depthM = depthM[:, :, 0:39]
    layerThicknessC = layerThicknessC[:, :, 0:39]
    layerThicknessM = layerThicknessM[:, :, 0:39]

#FW thickness (depth integrated)
VolFWC = np.nansum(np.multiply(layerThicknessC, FWC), axis=2)
VolFWM = np.nansum(np.multiply(layerThicknessM, FWM), axis=2)

#Index of maximum concentration
HmaxC = np.nanargmax(FWC[:,:,:], axis=2)
HmaxM = np.nanargmax(FWM[:,:,:], axis=2)

CmaxC = np.copy(lifwc[:,:,0])*0
DcmaxC = np.copy(lifwc[:,:,0])*0

CmaxM = np.copy(lifwc[:,:,0])*0
DcmaxM = np.copy(lifwc[:,:,0])*0

'''
Cmax_500 = np.copy(lifwc[:,:,0])*0
Dcmax_500 = np.copy(lifwc[:,:,0])*0

Cmax1000 = np.copy(lifwc[:,:,0])*0
Dcmax1000 = np.copy(lifwc[:,:,0])*0
'''

print('Looping')
for i in range(FWC.shape[0]):
    for j in range(FWC.shape[1]):
        CmaxC[i,j] = FWC[i,j,HmaxC[i,j]] #values of max concentration
        DcmaxC[i, j] = depthC[i, j, HmaxC[i, j]]  # depth of max concentration
        CmaxM[i, j] = FWM[i, j, HmaxM[i, j]]  # values of max concentration
        DcmaxM[i, j] = depthM[i, j, HmaxM[i, j]]  # depth of max concentration
        '''
        ind_500x = np.where(np.squeeze(depthC[i, j,:]) < 500)[0]
        if ind_500x.size > 0:
            HmaxC_500 = np.nanargmax(FWC[i, j, ind_500x])
            CmaxC[i, j] = FWC[i, j, HmaxC_500]
            DcmaxC[i, j] = depthC[i, j, HmaxC_500]
        else:
            CmaxC[i, j] = np.NaN
            DcmaxC[i, j] = np.NaN

        ind_500x = np.where(np.squeeze(depthM[i, j, :]) < 500)[0]
        if ind_500x.size > 0:
            HmaxC_500 = np.nanargmax(FWM[i, j, ind_500x])
            CmaxM[i, j] = FWM[i, j, HmaxC_500]
            DcmaxM[i, j] = depthM[i, j, HmaxC_500]
        else:
            CmaxM[i, j] = np.NaN
            DcmaxM[i, j] = np.NaN
        '''


print('DOne looping')


# Set up the figure and axis
ncols = 3
nrows = 3
fHeight = 12
fWidth = 18
cm = 1/2.54
fsize = 6

fig, ax = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.5*cm,wspace=0.1*cm)  # Increase vertical spacing

# Define custom colormap with 3 colors: white for 0, gray for NaN, transparent for 1
colors = [(1, 1, 1),  # White for 0
          (0.5, 0.5, 0.5),  # Gray for NaN
          (1, 1, 1, 0)]  # Transparent for 1
# Create a new colormap from the defined colors
binary_cmap = ListedColormap(colors, name='binary_cmap')
binary_cmap.set_over((1, 1, 1, 0))  # Make values < 1 (NaN) transparent
binary_cmap.set_bad((0.5, 0.5, 0.5))  # Set NaN to gray

cmap = plt.cm.RdYlBu_r

binary_mask = 0*DcmaxC+1
binary_mask[lat > -60] = 0

ctr = 0
for j in range(nrows):
    for k in range(ncols):
        if ctr == 0:
            ynow = VolFWC
            ttl = 'Freshwater thickness'
            cmin = 0
            cmax = 15
        elif ctr == 1:
            ynow = DcmaxC
            ttl = 'Depth of maximum concentration'
            cmin = 0
            if opt_nlay40 == 40:
                cmax = 1000
            else:
                cmax = 5000
        elif ctr == 2:
            ynow = CmaxC
            ttl = 'Maximum concentration'
            cmin = 0
            cmax = 0.02
        elif ctr == 3:
            ynow = VolFWM
            ttl = 'Freshwater thickness'
            cmin = 0
            cmax = 15
        elif ctr == 4:
            ynow = DcmaxM
            ttl = 'Depth of max concentration'
            cmin = 0
            if opt_nlay40 == 40:
                cmax = 1000
            else:
                cmax = 5000
        elif ctr == 5:
            ynow = CmaxM
            ttl = 'Max concentration'
            cmin = 0
            cmax = 0.02
        elif ctr == 6:
            ynow = VolFWM-VolFWC
            ttl = 'Freshwater thickness'
            cmap = plt.cm.seismic
            cmax = 3
            cmin = -cmax
        elif ctr == 7:
            ynow = DcmaxM-DcmaxC
            ttl = 'Depth of max concentration'
            cmax = 1000
            cmin = -cmax
        elif ctr == 8:
            ynow = CmaxM-CmaxC
            ttl = 'Max concentration'
            cmax = 0.005
            cmin = -cmax


        c = ax[j, k].imshow(ynow, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, alpha=1, vmin=cmin, vmax=cmax)
        ax[j,k].imshow(binary_mask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], cmap=binary_cmap, alpha=1)  # Adjust alpha for transparency

        #c = ax[j, k].pcolor(ynow, cmap=cmap, alpha=1)#, vmin=-cmax, vmax=cmax)
        #ax[j,k].pcolor(binary_mask, cmap=binary_cmap, alpha=1)  # Adjust alpha for transparency
        cbar = fig.colorbar(c, ax=ax[j, k], fraction=0.04, pad=0.04)
        cbar.ax.tick_params(labelsize=fsize - 1)  # Set tick label font size
        if j == 0:
            ax[j, k].set_title(ttl, fontsize=fsize-1)

        ax[j, k].autoscale(enable=True, axis='both', tight=True)
        ax[j, k].tick_params(axis='both', labelsize=fsize)
        ax[j, k].set_xlim([xmin, xmax])
        ax[j, k].set_ylim([ymin, ymax])
        ctr = ctr + 1
        if k == 0:
            ax[j, k].set_ylabel('Northing (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_yticklabels(), visible=False)
        if j == nrows - 1:
            ax[j, k].set_xlabel('Easting (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_xticklabels(), visible=False)

opt_save = 0
if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/FW_spatial/FW_spatial_diff_{fname}_{y1}_{y2}{nlayleg}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

fig = plt.figure()
DcmaxC[lat > -60] = np.NaN
DcmaxC1d = DcmaxC.reshape(-1)
DcmaxC1d = DcmaxC1d[~np.isnan(DcmaxC1d)]
plt.hist(DcmaxC1d, bins=1000, color='skyblue', edgecolor='black')
plt.show()