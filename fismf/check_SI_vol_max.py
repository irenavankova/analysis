#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

fpath = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)


si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/pismf/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'pismf_ave'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']
pismf_max_over_cells = siv.max(dim='nCells')
pismf_2_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=1)
pismf_10_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=9)

si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/fismf/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'fismf'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']

fismf_max_over_cells = siv.max(dim='nCells')
fismf_2_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=1)
fismf_10_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=9)

si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/hist/clim_2005-2014_ts_1951-2014/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'hist'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']

hist_max_over_cells = siv.max(dim='nCells')
hist_2_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=1)
hist_10_over_cells = siv.sortby('nCells', ascending=False).isel(nCells=9)



tseries_ds = xr.Dataset({
    'pismf_max_over_cells': pismf_max_over_cells,
    'fismf_max_over_cells': fismf_max_over_cells,
    'hist_max_over_cells': hist_max_over_cells,
    'pismf_2_over_cells': pismf_2_over_cells,
    'fismf_2_over_cells': fismf_2_over_cells,
    'hist_2_over_cells': hist_2_over_cells,
    'pismf_10_over_cells': pismf_10_over_cells,
    'fismf_10_over_cells': fismf_10_over_cells,
    'hist_10_over_cells': hist_10_over_cells
})

# Save to NetCDF
tseries_ds.to_netcdf(f'siv_max_over_cells.nc')

