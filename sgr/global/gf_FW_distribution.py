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

# MPAS Ocean mesh
fdir = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_101-110_ts_1-110/fw'
p_file = f'{fdir}/fw_9000.0x9000.0km_15.0km_Antarctic_stereo.nc'

#p_file = f"/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/remapped/layerThickness_6000.0x6000.0km_10.0km_Antarctic_stereo.nc"
#p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
ds = xr.open_dataset(p_file)
ssh = np.squeeze(ds['timeMonthly_avg_ssh'].values)
layerThickness = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)
lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
x = ds['x'].values
y = ds['y'].values

depth = layerThickness
depth[:,:,0] = depth[:,:,0] - ssh
depth = np.cumsum(depth, axis=2)
print(depth.shape)

plt.pcolor(np.squeeze(depth[:,:,0]))
plt.colorbar()
plt.show()

print(layerThickness.shape)
print([np.nanmin(layerThickness), np.nanmax(layerThickness)])

#plt.pcolor(np.squeeze(layerThickness[:,:,59]))
#plt.colorbar()
#plt.clim([0, 100])
#plt.show()

#conc = "/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_41-50_ts_1-50/conc_LIFW_41_50_6000.0x6000.0km_10.0km_Antarctic_stereo.nc"

# Load the NetCDF filec using xarray
#lf = xr.open_dataset(conc)
#lifwc = np.squeeze(lf['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
lifwc[lifwc < 0] = 0.0
lifwc = np.nan_to_num(lifwc, nan=0.0)
print(lifwc.shape)
print([np.nanmin(lifwc), np.nanmax(lifwc)])
print(np.where(lifwc > 1))
print(lifwc[np.where(lifwc < 0)])

#ind_x = np.where((x > 15) & (x < 19))
#plt.pcolor(np.squeeze(lifwc[:,:,1]))
#plt.colorbar()
#plt.show()
temp = np.multiply(layerThickness, lifwc)

print(temp.shape)
print([np.nanmin(temp), np.nanmax(temp)])

VolFW = np.nansum(temp, axis=2)

print(VolFW.shape)
print([np.nanmin(VolFW), np.nanmax(VolFW)])

'''
plt.pcolor(np.squeeze(VolFW))
plt.colorbar()
plt.show()
'''

HmaxC = np.nanargmax(lifwc[:,:,:], axis=2)

Cmax = np.copy(lifwc[:,:,0])*0
Dcmax = np.copy(lifwc[:,:,0])*0

print('Looping')
for i in range(lifwc.shape[0]):
    for j in range(lifwc.shape[1]):
        print(HmaxC[i,j])
        print(lifwc[i,j,HmaxC[i,j]])
        Cmax[i,j] = lifwc[i,j,HmaxC[i,j]] #values of max concentration
        print(depth[i, j, HmaxC[i, j]])
        Dcmax[i, j] = depth[i, j, HmaxC[i, j]]  # depth of max concentration

print('DOne looping')
plt.pcolor(np.squeeze(Dcmax))
plt.colorbar()
plt.show()

plt.pcolor(np.squeeze(Cmax))
plt.colorbar()
plt.show()

xx = 170
yy = 300

iz = HmaxC[xx, yy]
print(iz)
print(lifwc[xx, yy, iz])

plt.plot(np.arange(60),np.squeeze(lifwc[xx, yy, :]))
plt.plot(iz,lifwc[xx, yy, iz],'kx')
plt.show()


#HmaxC = np.argmax(lifwc, axis=2)

plt.pcolor(np.squeeze(HmaxC))
plt.colorbar()
plt.show()

'''
HmaxC = np.argmax(lifwc[:,:,30:59], axis=2)

plt.pcolor(np.squeeze(HmaxC))
plt.colorbar()
plt.show()
'''
