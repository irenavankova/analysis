#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

iceshelves = ["Antarctica", "Belli", "Amundsen", "Ross", "Eastant", "Amery", "Dml", "Fris", "Larsens"]
iam = gmask_reg.get_mask(iceshelves, mask_file)

# Step 2: Open the 100 NetCDF files and concatenate them along the 'time' dimension
ds = xr.open_mfdataset(f"{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc", combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True)
#ds = xr.open_mfdataset("data_*.nc", combine='by_coords', chunks={'time': 10})

print(ds.dims)

landIceFloatingMask = landIceFloatingMask.squeeze('Time')

# Melt rate
lifw = ds['timeMonthly_avg_landIceFreshwaterFlux'] 
lifw = lifw * landIceFloatingMask * areaCell * secPerYear / kgingt

# Step 4: Sum the masked temperature over the 'ncells' dimension for each time step
lifw_tseries_list = []

for n, shelf in enumerate(iceshelves):
    # Convert boolean mask to indices
    iis = np.where(iam[n])[0]

    # Apply index selection and sum
    tseries = lifw.isel(nCells=iis).sum(dim='nCells')
    lifw_tseries_list.append(tseries)

# Combine into a single DataArray with a new 'region' dimension
lifw_tseries = xr.concat(lifw_tseries_list, dim='region')
lifw_tseries['region'] = iceshelves  # label the new dimension
#for n in range(len(iceshelves)):
#    iis = iam[n, :]
#    lifw_tseries = lifw.isel(nCells=iis).sum(dim='nCells').compute()

print(lifw_tseries.dims)
lifw_tseries.name = 'lifw'

# Step 5: Save the resulting time series to a new NetCDF file

# Combine into a single Dataset
tseries_ds = xr.Dataset({
    'lifw': lifw_tseries,
})

# Save to NetCDF
tseries_ds.to_netcdf('lifw_reg_tseries.nc')

