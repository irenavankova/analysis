#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/pyremap/main/LICENSE

"""
Creates a mapping file that can be used with ncremap (NCO) to remap MPAS files
to a latitude/longitude grid.

Usage: Copy this script into the main MPAS-Analysis directory (up one level).
Modify the grid name, the path to the MPAS grid file and the output grid
resolution.
"""

import xarray

from pyremap import MpasCellMeshDescriptor, Remapper, get_polar_descriptor


# replace with the MPAS mesh name
inGridName = 'S123060'

# replace with the path to the desired mesh or restart file
# As an example, use:
# https://web.lcrc.anl.gov/public/e3sm/inputdata/ocn/mpas-o/oQU240/ocean.QU.240km.151209.nc

inFpath = '/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4'
inGridFileName = f'{inFpath}/ocean.SOwISC12to60E2r4.230220.nc'

inDescriptor = MpasCellMeshDescriptor(inGridFileName, inGridName)

# modify the size and resolution of the Antarctic grid as desired
outDescriptor = get_polar_descriptor(Lx=9000., Ly=9000., dx=15., dy=15.,
                                     projection='antarctic')
outGridName = outDescriptor.meshName

mappingFileName = f'{inFpath}/remap/map_SOwISC12to60E2r4_to_9000.0x9000.0km_15.0km_Antarctic_stereo_bilinear.nc'

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

# conservative remapping with 4 MPI tasks (using mpirun)
remapper.build_mapping_file(method='bilinear', mpiTasks=4)

# select the SST at the initial time as an example data set
#conc = "/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_41-50_ts_1-50/conc_LIFW_41_50.nc"
#p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
#fdir = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_101-110_ts_1-110/fw'
#conc = f'{fdir}/mpaso_ANN_010101_011012_climo.nc'
fdir = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_41-50_ts_1-50/fw'
conc = f'{fdir}/mpaso_ANN_004101_005012_climo.nc'


ds = xarray.open_dataset(conc)
dsOut = ds
#ds = xarray.open_dataset(p_file)
#dsOut = xarray.Dataset()
#dsOut['layerThickness'] = ds['layerThickness']

# do remapping again, this time with python remapping
print('remapping')
outFileName = f'{fdir}/fw_{outGridName}.nc'
#outFileName = f"/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_41-50_ts_1-50/conc_LIFW_41_50_{outGridName}.nc"
dsOut = remapper.remap(dsOut)
dsOut.to_netcdf(outFileName)
