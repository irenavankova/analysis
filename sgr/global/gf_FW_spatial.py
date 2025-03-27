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

# FW Ant projected file
fdir = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_101-110_ts_1-110/fw'
p_file = f'{fdir}/fw_9000.0x9000.0km_15.0km_Antarctic_stereo.nc'
ds = xr.open_dataset(p_file)

ssh = np.squeeze(ds['timeMonthly_avg_ssh'].values)
layerThickness = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)
lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
srfwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values)

mkm = 1000
x = ds['x'].values/mkm
y = ds['y'].values/mkm

depth = layerThickness
depth[:,:,0] = depth[:,:,0] - ssh
depth = np.cumsum(depth, axis=2)

lifwc[lifwc < 0] = 0.0
lifwc = np.nan_to_num(lifwc, nan=0.0)

srfwc[srfwc < 0] = 0.0
srfwc = np.nan_to_num(srfwc, nan=0.0)

#FW thickness (depth integrated)
VolFW = np.nansum(np.multiply(layerThickness, lifwc), axis=2)
VolFWsr = np.nansum(np.multiply(layerThickness, srfwc), axis=2)

#Index of maximum concentration
HmaxC = np.nanargmax(lifwc[:,:,:], axis=2)
HmaxCsr = np.nanargmax(srfwc[:,:,:], axis=2)


Cmax = np.copy(lifwc[:,:,0])*0
Dcmax = np.copy(lifwc[:,:,0])*0

Cmaxsr = np.copy(lifwc[:,:,0])*0
Dcmaxsr = np.copy(lifwc[:,:,0])*0

Cmax_500 = np.copy(lifwc[:,:,0])*0
Dcmax_500 = np.copy(lifwc[:,:,0])*0

Cmax1000 = np.copy(lifwc[:,:,0])*0
Dcmax1000 = np.copy(lifwc[:,:,0])*0

print('Looping')
for i in range(lifwc.shape[0]):
    for j in range(lifwc.shape[1]):
        Cmax[i,j] = lifwc[i,j,HmaxC[i,j]] #values of max concentration
        #if Cmax[i,j] >= 0.0001:
        Dcmax[i, j] = depth[i, j, HmaxC[i, j]]  # depth of max concentration
        #else:
        #    Dcmax[i, j] = np.NaN
        #ind_500 = np.where(np.squeeze(depth[i, j,:]) < 500)
        Cmaxsr[i, j] = srfwc[i, j, HmaxCsr[i, j]]  # values of max concentration
        # if Cmax[i,j] >= 0.0001:
        Dcmaxsr[i, j] = depth[i, j, HmaxCsr[i, j]]  # depth of max concentration

        '''
        ind_500x = np.where(np.squeeze(depth[i, j,:]) < 500)[0]
        #print(ind_500x)
        if ind_500x.size > 0:
            #ind_500 = ind_500x[0]

            #print(depth[i, j,ind_500])
            HmaxC_500 = np.nanargmax(lifwc[i, j, ind_500x])
            Cmax_500[i, j] = lifwc[i, j, HmaxC_500]
            if VolFW[i, j] >= 1:
                Dcmax_500[i, j] = depth[i, j, HmaxC_500]
            else:
                Dcmax_500[i, j] = np.NaN
        else:
            Cmax_500[i, j] = np.NaN
            Dcmax_500[i, j] = np.NaN

        '''

print('DOne looping')


# Set up the figure and axis
ncols = 3
nrows = 2
fHeight = 12
fWidth = 14
cm = 1/2.54

fig, ax = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.1*cm,wspace=0.1*cm)  # Increase vertical spacing

# Define custom colormap with 3 colors: white for 0, gray for NaN, transparent for 1
colors = [(1, 1, 1, 0),  # White for 0
          (0.5, 0.5, 0.5),  # Gray for NaN
          (1, 1, 1, 0)]  # Transparent for 1
# Create a new colormap from the defined colors
binary_cmap = ListedColormap(colors, name='binary_cmap')
binary_cmap.set_over((1, 1, 1, 0))  # Make values < 1 (NaN) transparent
binary_cmap.set_bad((0.5, 0.5, 0.5))  # Set NaN to gray

cmap = plt.cm.RdYlBu_r
fsize = 8

binary_mask = 0*Dcmax+1

ctr = 0
for j in range(nrows):
    for k in range(ncols):
        if ctr == 0:
            ynow = VolFW
            ttl = 'Freshwater thickness'
        elif ctr == 1:
            ynow = Dcmax
            ttl = 'Depth of maximum concentration'
        elif ctr == 2:
            ynow = Cmax
            ttl = 'Maximum concentration'
            cmax = 0.1
        elif ctr == 3:
            ynow = VolFWsr
            ttl = 'Freshwater thickness'
        elif ctr == 4:
            ynow = Dcmaxsr
            ttl = 'Depth of maximum concentration'
        elif ctr == 5:
            ynow = Cmaxsr
            ttl = 'Maximum concentration'
            cmax = 0.1

        c = ax[j, k].imshow(ynow, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, alpha=1)#, vmin=-cmax, vmax=cmax)
        ax[j,k].imshow(binary_mask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], cmap=binary_cmap, alpha=1)  # Adjust alpha for transparency

        #c = ax[j, k].pcolor(ynow, cmap=cmap, alpha=1)#, vmin=-cmax, vmax=cmax)
        #ax[j,k].pcolor(binary_mask, cmap=binary_cmap, alpha=1)  # Adjust alpha for transparency
        cbar = fig.colorbar(c, ax=ax[j, k], fraction=0.02, pad=0.04)
        cbar.ax.tick_params(labelsize=fsize - 1)  # Set tick label font size
        ax[j, k].set_title(ttl, fontsize=fsize-1)
        ax[j, k].autoscale(enable=True, axis='both', tight=True)
        ax[j, k].tick_params(axis='both', labelsize=fsize)
        #ax[j, k].set_xlim([xmin, xmax])
        #ax[j, k].set_ylim([ymin, ymax])
        ctr = ctr + 1
        if k == 0:
            ax[j, k].set_ylabel('Northing (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_yticklabels(), visible=False)
        if j == nrows - 1:
            ax[j, k].set_xlabel('Easting (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_xticklabels(), visible=False)


plt.show()