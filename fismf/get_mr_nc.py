#!/usr/bin/env python3

import xarray as xr
import numpy as np

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

# Step 2: Open the 100 NetCDF files and concatenate them along the 'time' dimension
ds = xr.open_mfdataset(f"{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc", combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True)
#ds = xr.open_mfdataset("data_*.nc", combine='by_coords', chunks={'time': 10})

print(ds.dims)

# Step 3: Apply the mask by multiplying the temperature variable with the mask
#landIceFloatingMask = landIceFloatingMask.expand_dims(time=ds['Time'], axis=0)
#areaCell = areaCell.expand_dims(time=ds['Time'], axis=0)
print(landIceFloatingMask.dims)
landIceFloatingMask = landIceFloatingMask.squeeze('Time')

print(areaCell.dims)
print(landIceFloatingMask.dims)
#landIceFloatingMask = landIceFloatingMask.squeeze('Time')
#print(landIceFloatingMask.dims)

# Melt rate
lifw = ds['timeMonthly_avg_landIceFreshwaterFlux'] 
lifw = lifw * landIceFloatingMask * areaCell
#lifw = lifw * areaCell

lifw = lifw * secPerYear / kgingt
print(lifw.dims)

# U* T*
Ti = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']
Tbl = ds['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature']
ustar = ds['timeMonthly_avg_landIceFrictionVelocity']

utstar = ustar * (Tbl - Ti)
utstar = utstar * landIceFloatingMask * areaCell
utstar = utstar * secPerYear / kgingt

# umag T

# Access variables
Vsq = ds['timeMonthly_avg_velocityMeridionalSquared']
#print(Vsq.dims)
Vsq = Vsq.isel(nVertLevels=0)
#print(Vsq.dims)
#Vsq = Vsq.squeeze()
#print(Vsq.dims)
#V = ds['timeMonthly_avg_velocityMeridional']
#V = V.isel(nVertLevels=0).squeeze()
Usq = ds['timeMonthly_avg_velocityZonalSquared']
#Usq = Usq.isel(nVertLevels=0).squeeze()
Usq = Usq.isel(nVertLevels=0)
#U = ds['timeMonthly_avg_velocityZonal']
#U = U.isel(nVertLevels=0).squeeze()
T = ds['timeMonthly_avg_activeTracers_temperature']
#T = T.isel(nVertLevels=0).squeeze()
T = T.isel(nVertLevels=0)

# Access the global attribute
Cd = ds.attrs['config_land_ice_flux_topDragCoeff']
utidal = ds.attrs['config_land_ice_flux_rms_tidal_velocity']

utsq = np.sqrt(Cd * (Usq + Vsq + utidal**2)) * (T - Ti)
#ut = np.sqrt(Cd * (U**2 + V**2 + utidal**2)) * (T - Ti)

utsq = utsq * landIceFloatingMask * areaCell
utsq = utsq * secPerYear / kgingt

#ut = ut * landIceFloatingMask * areaCell
#ut = ut * secPerYear / kgingt

# Step 4: Sum the masked temperature over the 'ncells' dimension for each time step
lifw_tseries = lifw.sum(dim='nCells').compute()
utstar_tseries = utstar.sum(dim='nCells').compute()
utsq_tseries = utsq.sum(dim='nCells').compute()

print(lifw_tseries.dims)
lifw_tseries.name = 'lifw'
print(utstar_tseries.dims)
utstar_tseries.name = 'utstar'
utsq_tseries.name = 'utsq'
# Step 5: Save the resulting time series to a new NetCDF file
#lifw_tseries.to_netcdf('lifw_tseries.nc')

# Combine into a single Dataset
tseries_ds = xr.Dataset({
    'lifw': lifw_tseries,
    'utstar': utstar_tseries,
    'utsq': utsq_tseries
})

# Save to NetCDF
tseries_ds.to_netcdf('proxy_tseries.nc')

