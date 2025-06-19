#!/usr/bin/env python3

import xarray as xr

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

# Step 2: Open the 100 NetCDF files and concatenate them along the 'time' dimension
ds = xr.open_mfdataset(f"{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.21*.nc", combine='by_coords')
#ds = xr.open_mfdataset("data_*.nc", combine='by_coords', chunks={'time': 10})

# Step 3: Apply the mask by multiplying the temperature variable with the mask
landIceFloatingMask = landIceFloatingMask.expand_dims(time=ds['Time'], axis=0)
areaCell = areaCell.expand_dims(time=ds['Time'], axis=0)

lifw = ds['timeMonthly_avg_landIceFreshwaterFlux'] * landIceFloatingMask * areaCell

lifw = lifw * secPerYear / kgingt

# Step 4: Sum the masked temperature over the 'ncells' dimension for each time step
lifw_tseries = lifw.sum(dim='ncells')

# Step 5: Save the resulting time series to a new NetCDF file
lifw_tseries.to_netcdf('lifw_tseries.nc')

