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
lnum = 7
#cmap_rho = 'gist_ncar'
cmap_conc_anom = 'cmo.balance'
cmap_conc = 'cmo.amp'

#cmap_rho = 'cmo.thermal'
vmax = 1
vmin = -vmax

plot_folder = f'{p_base}/plots'

fdir = f'{p_base}/sgr_on'
fdir_ref = f'{p_base}/sgr_off'

dsMesh = xarray.open_dataset(f'{fdir}/init.nc')
ds_ref = xarray.open_dataset(f'{fdir_ref}/timeSeriesStatsMonthly.0001-01-01.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0001-01-01.nc')

section_y = float(57000)
experiment = 'Ocean0'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

var_plot = 'lifw_vol'

ds.load()
lifw = np.squeeze(ds.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration.data)
thick = np.squeeze(ds.timeMonthly_avg_layerThickness.data)
lifw_thick_tot = np.nansum(lifw*thick, axis=1)
lifw_thick_tot = np.expand_dims(lifw_thick_tot, 0)
lifw_thick_all = lifw*thick
lifw_thick_L1 = lifw_thick_all[:, lnum]
lifw_thick_L1 = np.expand_dims(lifw_thick_L1, 0)

ds.timeMonthly_avg_landIceFreshwaterFlux.data = lifw_thick_tot

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_on', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

ds_ref.load()
lifw_ref = np.squeeze(ds_ref.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration.data)
thick_ref = np.squeeze(ds_ref.timeMonthly_avg_layerThickness.data)
lifw_thick_tot_ref = np.nansum(lifw_ref*thick_ref, axis=1)
lifw_thick_tot_ref = np.expand_dims(lifw_thick_tot_ref, 0)
lifw_thick_all_ref = lifw_ref*thick_ref
lifw_thick_L1_ref = lifw_thick_all_ref[:, lnum]
lifw_thick_L1_ref = np.expand_dims(lifw_thick_L1_ref, 0)

ds.timeMonthly_avg_landIceFreshwaterFlux.data = lifw_thick_tot_ref

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_off', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

ds.timeMonthly_avg_landIceFreshwaterFlux.data = lifw_thick_tot-lifw_thick_tot_ref

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_diff', oceanDomain='True',
    vmin=vmin/2, vmax=vmax/2, cmap=cmap_conc_anom,
    cmap_set_under='k', cmap_scale='linear')

var_plot = 'srfw_vol'

ds.load()
srfw = np.squeeze(ds.timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration.data)
thick = np.squeeze(ds.timeMonthly_avg_layerThickness.data)
srfw_thick_tot = np.nansum(srfw*thick, axis=1)
srfw_thick_tot = np.expand_dims(srfw_thick_tot, 0)

ds.timeMonthly_avg_landIceFreshwaterFlux.data = srfw_thick_tot

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_on', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

var_plot = f'lifw_volL{lnum}'
vmax = 0.5/2
vmin = -vmax

ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data*0
ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data + lifw_thick_L1

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_on', prefix=f'{var_plot}_on', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data*0
ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data + lifw_thick_L1_ref

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_off', prefix=f'{var_plot}_off', oceanDomain='True',
    vmin=0, vmax=vmax, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

ds.timeMonthly_avg_landIceFreshwaterFlux.data = lifw_thick_L1-lifw_thick_L1_ref

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle=f'{var_plot}_diff', prefix=f'{var_plot}_diff', oceanDomain='True',
    vmin=vmin/2, vmax=vmax/2, cmap=cmap_conc_anom,
    cmap_set_under='k', cmap_scale='linear')



'''
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
    
'''


