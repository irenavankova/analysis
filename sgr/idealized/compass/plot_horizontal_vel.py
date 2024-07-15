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
#cmap_rho = 'cmo.curl'
#cmap_rho = 'cmo.thermal'
cmap_rho = 'cmo.speed'
vmax = 3
vmin = 0

temp = 'rdn'
sgr = ["N", "A", "A"]
hloc = ["112", "132", "142"]

c = 0
out_name = f'{temp}_{hloc[c]}{sgr[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/horizontal/fvel/{cmap_rho}/{out_name}'

fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'

dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')

section_y = float(57000)
experiment = 'Ocean0'
var_plot = 'fvel'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

ds.timeMonthly_avg_landIceFrictionVelocity.data = ds.timeMonthly_avg_landIceFrictionVelocity.data * is_mask * 100

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFrictionVelocity,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}', oceanDomain='True',
    vmin=vmin, vmax=vmax, cmap=cmap_rho,
    cmap_set_under='k', cmap_scale='linear')


