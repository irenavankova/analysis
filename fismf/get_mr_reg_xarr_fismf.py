#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

fpath_pismf = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'
fpath_fismf = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370-fismf_ensmean/run/'

secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath_pismf}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file

mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

iceshelves = ["Antarctica", "Belli", "Amundsen", "Ross", "Eastant", "Amery", "Dml", "Fris", "Larsens"]
iam = gmask_reg.get_mask(iceshelves, mask_file)

# Step 2: Open the 100 NetCDF files and concatenate them along the 'time' dimension
ds = xr.open_mfdataset(f"{fpath_fismf}v2_1.SORRM.ssp370-fismf_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc", combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True)

landIceFloatingMask = landIceFloatingMask.squeeze('Time')

# Melt rate
lifw = ds['timeMonthly_avg_landIceFreshwaterFlux'] 
lifw = lifw * landIceFloatingMask * areaCell * secPerYear / kgingt

# --- Vectorized masking over all regions ---
# Stack region mask into xarray DataArray for better alignment
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': iceshelves})

# Expand lifw to shape (region, Time, nCells) by broadcasting
lifw_broadcasted = lifw.expand_dims(region=iam_da.region)

# Use masking across all regions simultaneously
lifw_masked = lifw_broadcasted.where(iam_da)

# Sum over nCells → result shape: (region, Time)
lifw_tseries = lifw_masked.sum(dim='nCells')

# Transpose to have (Time, region) if desired
lifw_tseries = lifw_tseries.transpose('Time', 'region')
lifw_tseries.name = 'lifw'

# Save to NetCDF
lifw_tseries.to_dataset().to_netcdf('lifw_reg_fismf_ave_tseries.nc')

