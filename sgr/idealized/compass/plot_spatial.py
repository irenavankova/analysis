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

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'

temp = 'rdn'
sgr = ["R", "A", "A", "A", "A", "B", "B"]
hloc = ["112", "112", "132", "122", "142", "132", "142"]

c = 4
out_name = f'{temp}_{hloc[c]}{sgr[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/spatial/{out_name}'

fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'

dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')

section_y = float(40000)

experiment = 'Ocean0'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

vmax = 80
rho_fw = 1000
sec_per_year = 86400. * 365.
ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data / rho_fw * sec_per_year
ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data * is_mask

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFreshwaterFlux,
    nameInTitle='landIceFreshwaterFlux', prefix='melt', oceanDomain='True',
    vmin=-vmax, vmax=vmax, cmap='cmo.balance',
    cmap_set_under='k', cmap_scale='linear')


maxLevelCell = dsMesh.maxLevelCell - 1
dsBot = xarray.DataArray(ds.timeMonthly_avg_potentialDensity)
dsBot.coords['verticalIndex'] = ('nVertLevels', np.arange(dsBot.sizes['nVertLevels']))
dsBot = dsBot.where(dsBot.verticalIndex == maxLevelCell)
dsBot = dsBot.sum(dim='nVertLevels').where(maxLevelCell >= 0)

plotter.plot_horiz_series(
    dsBot,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_bot', oceanDomain='True',
    vmin=1027.0, vmax=1028.0,
    cmap_set_under='k', cmap_scale='linear')

minLevelCell = dsMesh.minLevelCell - 1
dsTop = xarray.DataArray(ds.timeMonthly_avg_potentialDensity)
dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsTop,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_top', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear')

#dsDiff = xarray.DataArray(dsTop)
dsDiff = dsTop-dsBot

plotter.plot_horiz_series(
    dsDiff,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_diff', oceanDomain='True',
    vmin=-0.8, vmax=0,
    cmap_set_under='k', cmap_scale='linear')

section_y = float(25000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_vert25', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap='gist_ncar')

section_y = float(40000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_vert40', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap='gist_ncar')

section_y = float(55000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_vert55', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap='gist_ncar')