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

y1 = 101
y2 = 110

ncname = f'mpaso_ANN_{str(y1).zfill(4)}01_{str(y2).zfill(4)}12_climo.nc'
fdirC = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{y1}-{y2}_ts_1-{y2}/anytr/netcdf_SOwISC12to60E2r4_to_up-Larsen-C/{ncname}'
fdirM = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{y1}-{y2}_ts_1-{y2}/anytr/netcdf_SOwISC12to60E2r4_to_up-Larsen-C/{ncname}'

# Load the NetCDF filec using xarray
dC = xr.open_dataset(fdirC)
dM = xr.open_dataset(fdirM)

lifw = np.squeeze(dM['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values) - np.squeeze(dC['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
sifw = np.squeeze(dM['timeMonthly_avg_freshwaterTracers_seaIceFreshWaterConcentration'].values) - np.squeeze(dC['timeMonthly_avg_freshwaterTracers_seaIceFreshWaterConcentration'].values)
srfw = np.squeeze(dM['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values) - np.squeeze(dC['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values)

#x = lf['dNode'].values
#z = lf['y'].values
cmap = plt.cm.seismic
cmax = 0.0005

ncols = 2
nrows = 2
fHeight = 8
fWidth = 16
cm = 1/2.54
fig, ax = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
ax[0,0].pcolor(np.flipud(lifw.T), cmap=cmap, vmin=-cmax, vmax=cmax)
ax[0,1].pcolor(np.flipud(sifw.T), cmap=cmap, vmin=-cmax, vmax=cmax)
ax[1,0].pcolor(np.flipud(srfw.T), cmap=cmap, vmin=-cmax, vmax=cmax)
ax[1,1].pcolor(np.flipud(lifw.T + sifw.T+ srfw.T), cmap=cmap, vmin=-cmax, vmax=cmax)
plt.show()


