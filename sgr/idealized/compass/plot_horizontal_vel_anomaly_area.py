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
cmap_rho = 'cmo.balance'
#cmap_rho = 'cmo.thermal'
vmax = 1
vmin = -vmax

d_scale = 30000

temp = 'rdn'
#sgr = ["A", "A","A", "A"]
#hloc = ["112", "132", "122", "142"]
sgr = ["A", "A", "B", "B", "D", "D", "C", "C", "R", "R"]
hloc = ["122", "142", "122", "142", "122", "142", "122", "142", "122", "142"]

c = 0
out_name = f'{temp}_{hloc[c]}{sgr[c]}'
plot_folder = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/horizontal/fvel_anom_area/{out_name}'

fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'
fdir_ref = f'{p_base}/{temp}/{temp}_112N'

dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
ds_ref = xarray.open_dataset(f'{fdir_ref}/timeSeriesStatsMonthly.0002-12-01.nc')
ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')

x = np.squeeze(dsMesh.xCell.data)
y = np.squeeze(dsMesh.yCell.data)
z = np.squeeze(dsMesh.zCell.data)

fdirf = f'{p_base}/rd/rd_{hloc[c]}R'
df = xarray.open_dataset(f'{fdirf}/sgr_data.nc')
df.load()
subglacialRunoffFlux = np.squeeze(df.subglacialRunoffFlux.data)
sgr_iii = subglacialRunoffFlux > 0
sgr_pt = np.where(sgr_iii)[0]

section_y = float(57000)
experiment = 'Ocean0'
var_plot = 'fvel_anom'
plotter = MoviePlotter(inFolder=work_dir,
                       streamfunctionFolder=work_dir,
                       outFolder=plot_folder, sectionY=section_y,
                       dsMesh=dsMesh, ds=ds, expt=experiment,
                       showProgress=False)

is_mask = dsMesh.landIceFloatingMask.data
is_mask = np.where(is_mask==0, np.nan, is_mask)

dist2point = np.absolute(np.sqrt((x - x[sgr_pt]) ** 2 + (y - y[sgr_pt]) ** 2 + (z - z[sgr_pt]) ** 2))
# ii = numpy.where(dist2point == dist2point.min())
i_ave = dist2point < d_scale
i_ave_mask = np.where(dist2point > d_scale, np.nan, is_mask)

ts = ds.timeMonthly_avg_landIceFrictionVelocity.data-ds_ref.timeMonthly_avg_landIceFrictionVelocity.data

ds.timeMonthly_avg_landIceFrictionVelocity.data = ts * is_mask * 100 * i_ave_mask

plotter.plot_horiz_series(
    ds.timeMonthly_avg_landIceFrictionVelocity,
    nameInTitle=f'{var_plot}', prefix=f'{var_plot}', oceanDomain='True',
    vmin=vmin, vmax=vmax, cmap=cmap_rho,
    cmap_set_under='k', cmap_scale='linear')


