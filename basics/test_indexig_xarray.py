#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xarray

d2y = 365

tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
#tseries = ["0001-0050", "0001-0050", "0001-0050", "0001-0050"]
#tsegment = ["clim_41-50_ts_1-50", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
sims = ["S12_control", "S12_mali", "S12_mali_x4", "S12_mali_x8"]
legsims = ["control", "sgr-x1", "sgr-x4" , "sgr-x8"]

s = 3
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/iceShelfFluxes_{tseries[s]}.nc'
dsOut = xarray.open_dataset(p_file)
dsOut.load()
#dsOut.isel(nRegions=1)
dsOut = dsOut.sel(nRegions=[2, 4],Time=dsOut.Time[dsOut.Time <= d2y*25])

regionNames = dsOut.regionNames.data
print(dsOut)
regionNames = np.squeeze(regionNames[0,:])
print(regionNames)
#Time = np.squeeze(dsOut.Time.data) / d2y









