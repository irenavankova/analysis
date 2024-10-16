#!/usr/bin/env python3
import os
import shutil
import time

import numpy as np
import xarray
import cmocean  # noqa: F401


from plot_iv import MoviePlotter

# get the current working directory
current_working_directory = os.getcwd()
# print output to the console
print(current_working_directory)
work_dir = current_working_directory

p_base = '/Users/irenavankova/Work/data_sim/SGR/test_fw_tracers/test_CB_FW_mpaso/'

#cmap_rho = 'gist_ncar'
cmap_conc_anom = 'cmo.balance'
cmap_conc = 'cmo.amp'

#cmap_rho = 'cmo.thermal'
vmax = 0.05
vmin = -vmax

plot_folder = f'{p_base}/plots'

fdir = f'{p_base}/sgr_on'
fdir_ref = f'{p_base}/sgr_off'

dsMesh = xarray.open_dataset(f'{fdir}/init.nc')
ds_ref = xarray.open_dataset(f'{fdir_ref}/timeSeriesStatsMonthly.0001-01-01.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0001-01-01.nc')

section_y = float(57000)
experiment = 'Ocean0'
var_plot = 'lifw'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

maxLevelCell = dsMesh.maxLevelCell - 1
minLevelCell = dsMesh.minLevelCell - 1

dsTop = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration)
dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsTop,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_on', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

dsRefTop = xarray.DataArray(ds_ref.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration)
dsRefTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsRefTop.sizes['nVertLevels']))
dsRefTop = dsRefTop.where(dsRefTop.verticalIndex == minLevelCell)
dsRefTop = dsRefTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsRefTop,
    nameInTitle=f'{var_plot}_off', prefix=f'{var_plot}_off', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')


dsDiffTop = dsTop - dsRefTop
plotter.plot_horiz_series(
    dsDiffTop,
    nameInTitle=f'{var_plot}_diff', prefix=f'{var_plot}_diff', oceanDomain='True',
    vmin=vmin/2, vmax=vmax/2, cmap=cmap_conc_anom,
    cmap_set_under='k', cmap_scale='linear')

var_plot = 'srfw'

dsTop = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration)
dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsTop,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_on', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')



