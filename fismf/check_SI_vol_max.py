#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg


si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/pismf/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'pismf_ave'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']
pismf_max_over_cells = siv.max(dim='nCells')

si_fname = '/lcrc/group/e3sm/ac.vankova/scratch/mpas_analysis/FISMF/fismf_701/clim_2091-2100_ts_2015-2100/timeseries/mpasTimeSeriesSeaIce.nc'
out_fname = 'fismf_701'
ds = xr.open_dataset(si_fname, chunks={'time': 12}, decode_timedelta=True)

siv = ds['timeMonthly_avg_iceVolumeCell']
fismf_max_over_cells = siv.max(dim='nCells')


tseries_ds = xr.Dataset({
    'pismf_max_over_cells': pismf_max_over_cells,
    'fismf_max_over_cells': fismf_max_over_cells
})

# Save to NetCDF
tseries_ds.to_netcdf(f'siv_max_over_cells.nc')

