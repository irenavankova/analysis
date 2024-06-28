#!/usr/bin/env python3
import os
import shutil
import time

import numpy as np
import xarray

from compass.model import run_model
from compass.ocean.tests.isomip_plus.evap import update_evaporation_flux
from compass.ocean.tests.isomip_plus.viz.plot import MoviePlotter
from compass.step import Step

# seconds per day
sec_per_day = 86400.
# seconds in a year
sec_per_year = sec_per_day * 365.
# seconds per month (approximate)
sec_per_month = sec_per_day * 365. / 12.
# density of freshwater (kg/m^3)
rho_fw = 1000.

# get the current working directory
current_working_directory = os.getcwd()
# print output to the console
print(current_working_directory)
work_dir = current_working_directory

#Get location of data
#case_num = ["ra", "rb", "rc", "rd", "re", "rg"]
case_num = ["ra"]
Tbot = np.array([-1.9, -1, 0, 1, 2, 4])
loc_sgr = 112
#loc_case = ["N", "R", "A", "C", "B"]
loc_case = ["A"]
sgr = np.array([0, 1, 10, 50, 100])
n_case = len(case_num)
n_sgr = len(sgr)

for c in range(n_case):
    for s in range(n_sgr):

        work_dir = f'/lcrc/group/e3sm/ac.vankova/compass/sg_tests/sg_pull_w_fraz_yesC/{case_num[c]}/{case_num[c]}_{loc_sgr}{loc_case[s]}/ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation'

        mesh_array = f'{work_dir}/init.nc'
        experiment = 'Ocean0'
        dsMesh = xarray.open_dataset(mesh_array)

#for yr in range(1, 3, 1):
#    yy = f'{yr}'.zfill(2)
#    for k in range(1, 13, 1):
        #Get snapshots from restart files
        yy = f'{2}'.zfill(2)
        mm = f'{12}'.zfill(2)
        output_array = f'{work_dir}/timeSeriesStatsMonthly.00{yy}-{mm}-01.nc'
        if os.path.exists(output_array):
            print(output_array)

            # plot a few fields
            plot_folder = f'{work_dir}/plots/final'

            ds = xarray.open_dataset(output_array)

            section_y = float(40000)
            # show progress only if we're not writing to a log file
            plotter = MoviePlotter(inFolder=work_dir,
                                   streamfunctionFolder=work_dir,
                                   outFolder=plot_folder, sectionY=section_y,
                                   dsMesh=dsMesh, ds=ds, expt=experiment,
                                   showProgress=False)
            '''
            plotter.plot_3d_field_top_bot_section(
                ds.timeMonthly_avg_potentialDensity, nameInTitle='potential density', prefix='pdens', suffix=f'{yy}-{mm}',
                units='kg/m^3', vmin=1027.2, vmax=1027.8, cmap='cmo.thermal')

            plotter.plot_3d_field_top_bot_section(
                ds.timeMonthly_avg_activeTracers_temperature, nameInTitle='temperature', prefix='temp', suffix=f'{yy}-{mm}',
                units='C', vmin=-2., vmax=1., cmap='cmo.thermal')

            plotter.plot_3d_field_top_bot_section(
                ds.timeMonthly_avg_activeTracers_salinity, nameInTitle='salinity', prefix='salin', suffix=f'{yy}-{mm}',
                units='PSU', vmin=33.8, vmax=34.7, cmap='cmo.thermal')

            if 'timeMonthly_avg_kineticEnergyCell' in ds.keys():
                plotter.plot_3d_field_top_bot_section(
                    ds.timeMonthly_avg_kineticEnergyCell, nameInTitle='KECell', prefix='KE', suffix=f'{yy}-{mm}',
                    units='PSU', vmin=33.8, vmax=34.7, cmap='cmo.thermal')
            '''

            # constants.sec_per_year / constants.rho_fw

            if 'timeMonthly_avg_landIceFreshwaterFlux' in ds.keys():
                plotter.plot_horiz_series(
                    ds.timeMonthly_avg_landIceFreshwaterFlux,
                    nameInTitle='landIceFreshwaterFlux', prefix='liff', suffix=f'{yy}-{mm}', oceanDomain='True',
                    vmin=0, vmax=1e-4,
                    cmap_set_under='k', cmap_scale='linear')



            #plotter.plot_3d_field_top_bot_section(
            #    ds.timeMonthly_avg_activeTracers_salinity, nameInTitle='salinity', prefix='salin', suffix=f'{yy}-{mm}',
            #    units='PSU', vmin=33.8, vmax=34.7, cmap='cmo.haline')

            #plotter.plot_3d_field_top_bot_section(
            #    ds.timeMonthly_avg_vertVelocityTop, nameInTitle='vertical velocity', prefix='vertVel', suffix=f'{yy}-{mm}',
            #    units='m/s', vmin=-0.000001, vmax=0.000001, cmap='cmo.haline')


