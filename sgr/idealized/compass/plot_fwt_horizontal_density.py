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

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sgr_tracers_mpaso'

#cmap_rho = 'gist_ncar'
#cmap_rho = 'cmo.dense'
cmap_rho = 'Spectral_r'
cmap_conc = 'cmo.amp'

temp = 'FS'
hloc = ["011", "001", "010", "110", "111"]

c = 1
out_name = f'{temp}_{hloc[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/fwt/horizontal/{out_name}'

fdir = f'{p_base}/{temp}_{hloc[c]}'

dsMesh = xarray.open_dataset(f'{fdir}/restart.0002-01-01_00.00.00.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0001-12-01.nc')

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

suff = 55
section_y = float(suff*1000)
experiment = 'Ocean0'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

try:
    vmax = 80.0
    rho_fw = 1000
    sec_per_year = 86400. * 365.
    ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data / rho_fw * sec_per_year
    ds.timeMonthly_avg_landIceFreshwaterFlux.data = ds.timeMonthly_avg_landIceFreshwaterFlux.data * is_mask

    plotter.plot_horiz_series(
        ds.timeMonthly_avg_landIceFreshwaterFlux,
        nameInTitle='landIceFreshwaterFlux', prefix='melt', oceanDomain='True',
        vmin=-vmax, vmax=vmax, cmap='cmo.balance',
        cmap_set_under='k', cmap_scale='linear')
except:
    print("Melt not activated")

maxLevelCell = dsMesh.maxLevelCell - 1
minLevelCell = dsMesh.minLevelCell - 1

dsBot = xarray.DataArray(ds.timeMonthly_avg_potentialDensity)
dsBot.coords['verticalIndex'] = ('nVertLevels', np.arange(dsBot.sizes['nVertLevels']))
dsBot = dsBot.where(dsBot.verticalIndex == maxLevelCell)
dsBot = dsBot.sum(dim='nVertLevels').where(maxLevelCell >= 0)

plotter.plot_horiz_series(
    dsBot,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_bot', oceanDomain='True',
    vmin=1027.0, vmax=1030.0, cmap=cmap_rho,
    cmap_set_under='k', cmap_scale='linear')

dsTop = xarray.DataArray(ds.timeMonthly_avg_potentialDensity)
dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsTop,
    nameInTitle='landIceFreshwaterFlux', prefix='rho_top', oceanDomain='True',
    vmin=1027, vmax=1030, cmap=cmap_rho,
    cmap_set_under='k', cmap_scale='linear')

dsBot = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration)
dsBot.coords['verticalIndex'] = ('nVertLevels', np.arange(dsBot.sizes['nVertLevels']))
dsBot = dsBot.where(dsBot.verticalIndex == maxLevelCell)
dsBot = dsBot.sum(dim='nVertLevels').where(maxLevelCell >= 0)

plotter.plot_horiz_series(
    dsBot,
    nameInTitle='landIceFreshwaterFlux', prefix='csg_bot', oceanDomain='True',
    vmin=0, vmax=0.1, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

dsTop = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration)
dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

plotter.plot_horiz_series(
    dsTop,
    nameInTitle='landIceFreshwaterFlux', prefix='csg_top', oceanDomain='True',
    vmin=0, vmax=0.1, cmap=cmap_conc,
    cmap_set_under='k', cmap_scale='linear')

try:
    dsBot = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration)
    dsBot.coords['verticalIndex'] = ('nVertLevels', np.arange(dsBot.sizes['nVertLevels']))
    dsBot = dsBot.where(dsBot.verticalIndex == maxLevelCell)
    dsBot = dsBot.sum(dim='nVertLevels').where(maxLevelCell >= 0)

    plotter.plot_horiz_series(
        dsBot,
        nameInTitle='landIceFreshwaterFlux', prefix='cli_bot', oceanDomain='True',
        vmin=0, vmax=0.1, cmap=cmap_conc,
        cmap_set_under='k', cmap_scale='linear')

    dsTop = xarray.DataArray(ds.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration)
    dsTop.coords['verticalIndex'] = ('nVertLevels', np.arange(dsTop.sizes['nVertLevels']))
    dsTop = dsTop.where(dsTop.verticalIndex == minLevelCell)
    dsTop = dsTop.sum(dim='nVertLevels').where(minLevelCell >= 0)

    plotter.plot_horiz_series(
        dsTop,
        nameInTitle='landIceFreshwaterFlux', prefix='cli_top', oceanDomain='True',
        vmin=0, vmax=0.1, cmap=cmap_conc,
        cmap_set_under='k', cmap_scale='linear')
except:
    print("Melt tracer not activated")


'''

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='rho', prefix=f'rho_vert{suff}', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration,
    nameInTitle='csg', prefix=f'csg_vert{suff}', oceanDomain='True',
    vmin=0, vmax=0.1,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_conc)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration,
    nameInTitle='cli', prefix=f'cli_vert{suff}', oceanDomain='True',
    vmin=0, vmax=0.1,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_conc)

section_y = float(25000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='rho', prefix='rho_vert25', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

section_y = float(55000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='rho', prefix='rho_vert55', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)
'''