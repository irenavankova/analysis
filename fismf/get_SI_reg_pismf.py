#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/pismf/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

iceshelves = ["si_ab", "si_rossw", "si_rosse", "si_so", "si_ea", "si_amery", "si_dml", "si_weddell"]

iam = gmask_reg.get_mask(iceshelves, mask_file)

ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

landIceFloatingMask = landIceFloatingMask.squeeze('Time')

sic = ds['timeMonthly_avg_iceAreaCell']

# --- Vectorized masking over all regions ---
# Stack region mask into xarray DataArray for better alignment
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': iceshelves})

# Expand to shape (region, Time, nCells) by broadcasting
sic_broadcasted = sic.expand_dims(region=iam_da.region)
areaCell_broadcasted = areaCell.expand_dims(region=iam_da.region)

# Use masking across all regions simultaneously
sic_masked = sic_broadcasted.where(iam_da)
areaCell_masked = areaCell_broadcasted.where(iam_da)

sic_weighted = sic_masked * areaCell_masked

# Sum over nCells → result shape: (region, Time)
sic_tseries = sic_weighted.sum(dim='nCells')/areaCell_masked.sum(dim='nCells')

# Transpose to have (Time, region) if desired
sic_tseries = sic_tseries.transpose('Time', 'region')
sic_tseries.name = 'sic'

# Save to NetCDF
sic_tseries.to_dataset().to_netcdf('sic_reg_pismf_ave_tseries.nc')

