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
cmap_rho = 'Spectral_r'

temp = 'rd'
sgr = ["N", "A", "A"]
hloc = ["112", "132", "142"]

c = 0
out_name = f'{temp}_{hloc[c]}{sgr[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/vertical/{cmap_rho}/{out_name}'

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

plotter.plot_vertical_section(
    ds.timeMonthly_avg_potentialDensity,
    nameInTitle='rho', prefix='rho_vert40', oceanDomain='True',
    vmin=1027, vmax=1028,
    cmap_set_under='k', cmap_scale='linear', cmap=cmap_rho)

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