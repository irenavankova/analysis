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
siv = siv/areaCell
pismf_max_over_cells = siv.max(dim='nCells')

si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/fismf_701/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'fismf_701'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']
siv = siv/areaCell

fismf_max_over_cells = siv.max(dim='nCells')


tseries_ds = xr.Dataset({
    'pismf_max_over_cells': pismf_max_over_cells,
    'fismf_max_over_cells': fismf_max_over_cells
})

# Save to NetCDF
tseries_ds.to_netcdf(f'siv_max_over_cells.nc')

