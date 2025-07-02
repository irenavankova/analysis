#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

out_fname = 'pismf_ave'
#out_fname = 'hist_ave'
#out_fname = 'fismf_701'

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

if out_fname == 'pismf_ave':
    fnc = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/v2_1.SORRM.ssp370_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc'
elif out_fname == 'hist_ave':
    fnc = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.historical_ensmean/run/v2_1.SORRM.historical_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc'
elif out_fname == 'fismf_701':
    ymem = 701
    fpath_fismf = f'/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_0{ymem}-fismf/'
    if ymem == 701:
        fpath_fismf = f'{fpath_fismf}archive/ocn/hist/'
    elif ymem == 751:
        fpath_fismf = f'{fpath_fismf}run/'
    fnc = f"{fpath_fismf}v2_1.SORRM.ssp370_0{ymem}-fismf.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

iceshelves = ["si_so", "si_ab", "si_rosse", "si_rossw", "si_ea", "si_amery", "si_dml", "si_weddell"]

iam = gmask_reg.get_mask(iceshelves, mask_file)

ds = xr.open_mfdataset(fnc, combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True)

tau_zonal = ds['timeMonthly_avg_windStressZonal']
tau_merid = ds['timeMonthly_avg_windStressMeridional']

# --- Vectorized masking over all regions ---
# Stack region mask into xarray DataArray for better alignment
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': iceshelves})

# Expand to shape (region, Time, nCells) by broadcasting
tau_zonal_broadcasted = tau_zonal.expand_dims(region=iam_da.region)
tau_merid_broadcasted = tau_merid.expand_dims(region=iam_da.region)

areaCell_broadcasted = areaCell.expand_dims(region=iam_da.region)

# Use masking across all regions simultaneously
tau_zonal_masked = tau_zonal_broadcasted.where(iam_da)
tau_merid_masked = tau_merid_broadcasted.where(iam_da)

areaCell_masked = areaCell_broadcasted.where(iam_da)

tau_zonal_weighted = tau_zonal_masked * areaCell_masked
tau_merid_weighted = tau_merid_masked * areaCell_masked

# Sum over nCells → result shape: (region, Time)
tau_zonal_tseries = tau_zonal_weighted.sum(dim='nCells')/areaCell_masked.sum(dim='nCells')
tau_merid_tseries = tau_merid_weighted.sum(dim='nCells')/areaCell_masked.sum(dim='nCells')

# Transpose to have (Time, region) if desired
tau_zonal_tseries = tau_zonal_tseries.transpose('Time', 'region')
tau_zonal_tseries.name = 'tau_zonal'
tau_merid_tseries = tau_merid_tseries.transpose('Time', 'region')
tau_merid_tseries.name = 'tau_merid'

tseries_ds = xr.Dataset({
    'tau_zonal': tau_zonal_tseries,
    'tau_merid': tau_merid_tseries
})

# Save to NetCDF
tseries_ds.to_netcdf(f'tau_reg_{out_fname}_tseries.nc')

