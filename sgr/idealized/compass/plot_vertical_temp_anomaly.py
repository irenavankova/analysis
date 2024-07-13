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

#cmap_rho = 'gist_ncar'
#cmap_rho = 'cmo.dense'
cmap_rho = 'cmo.balance'
vmin = -0.2
vmax = -vmin

temp = 'rdn'
sgr = ["A", "A", "A", "A"]
hloc = ["112", "132", "122", "142"]

c = 0
out_name = f'{temp}_{hloc[c]}{sgr[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/vertical/temp_anom/{cmap_rho}/{out_name}'

fdir_ref = f'{p_base}/{temp}/{temp}_112N'
fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'

dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
ds_ref = xarray.open_dataset(f'{fdir_ref}/timeSeriesStatsMonthly.0002-12-01.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')

tref = ds_ref.timeMonthly_avg_activeTracers_temperature.data
ds.timeMonthly_avg_activeTracers_temperature.data = ds.timeMonthly_avg_activeTracers_temperature.data-tref

#ds_anom = ds-ds_ref

section_y = float(57000)
experiment = 'Ocean0'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)
var_plot = 'temp'
plotter.plot_vertical_section(
    ds.timeMonthly_avg_activeTracers_temperature,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}_vert57', oceanDomain='True',
    vmin=vmin, vmax=vmax,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

section_y = float(40000)
experiment = 'Ocean0'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)
var_plot = 'temp'
plotter.plot_vertical_section(
    ds.timeMonthly_avg_activeTracers_temperature,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}_vert40', oceanDomain='True',
    vmin=vmin, vmax=vmax,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

section_y = float(25000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_activeTracers_temperature,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}_vert25', oceanDomain='True',
    vmin=vmin, vmax=vmax,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

section_y = float(55000)
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

plotter.plot_vertical_section(
    ds.timeMonthly_avg_activeTracers_temperature,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}_vert55', oceanDomain='True',
    vmin=vmin, vmax=vmax,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)