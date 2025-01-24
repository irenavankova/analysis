#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import gmask_is

# Get subglacial discharge per ice shelf
rho_fw = 1000
d2y = 365

#iceshelves = ["Amery", "Thwaites", "Getz", "Totten"]
iceshelves = ["George_VI", "Pine_Island", "Thwaites", "Getz", "Ross", "Totten", "Amery","Fimbul","Filchner-Ronne", "Larsen_C"]

sgr_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/sgr_files_IV/MALI_20240327/DSGR.massFlux.MALI.out2055.SOwISC12to60E2r4.20240328.nc'
ds = xarray.open_dataset(sgr_file)
ds.load()
sgr = np.squeeze(ds.subglacialRunoffFlux.data)

iam, areaCell = gmask_is.get_mask(iceshelves)

sgr_vol_flux = np.zeros((len(iceshelves)))
for n in range(len(iceshelves)):
    iis = iam[n, :]
    sgr_mass_flux = np.nansum(sgr[iis] * areaCell[iis], axis=0)
    sgr_vol_flux[n] = sgr_mass_flux/rho_fw
    print(sgr_vol_flux[n])

# Get melt rates
tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
sims = ["S12_control", "S12_mali", "S12_mali_x4" , "S12_mali_x8"]

melt = np.zeros((len(iceshelves), len(sims)))

for n in range(len(iceshelves)):
    shelf_name = iceshelves[n]
    if shelf_name == "PineIsland":
        shelf_name = "Pine_Island"

    print(shelf_name)
    for s in range(len(sims)):

        p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/iceShelfFluxes_{tseries[s]}.nc'

        dsOut = xarray.open_dataset(p_file)
        dsOut.load()
        MeltFluxNow = np.squeeze(dsOut.integratedMeltFlux.data)
        Time = np.squeeze(dsOut.Time.data) / d2y
        if s == 0:
            regionNames = dsOut.regionNames.data
            regionNames = np.squeeze(regionNames[0, :])
            ind_x = np.where((regionNames == shelf_name))
            print(regionNames[ind_x])

        ind_t = np.where((Time >= 21) & (Time <= 26))
        melt[n,s] = np.mean(MeltFluxNow[ind_t, ind_x])

fHeight = 5
fWidth = 10
plt.figure(figsize=(fWidth, fHeight))

for n in range(len(iceshelves)):
    #sgr_xax = np.array([0, 1, 4, 8]) * sgr_vol_flux[n]
    #plt.plot(sgr_xax, melt[n,:], linewidth=1, label = iceshelves[n])
    sgr_xax = np.array([0, 1, 4, 8]) * sgr_vol_flux[n]
    plt.plot(sgr_xax, melt[n,:]-melt[n,0], 'x', linewidth=1, linestyle = '-', label = iceshelves[n])

plt.xlabel('SGR (m^3/s)')
plt.ylabel('Integrated melf flux (GT/a)')
plt.legend(loc = 2,fontsize=8)
plt.grid()
plt.show()

fHeight = 5
fWidth = 10
plt.figure(figsize=(fWidth, fHeight))

for n in range(len(iceshelves)):
    #sgr_xax = np.array([0, 1, 4, 8]) * sgr_vol_flux[n]
    #plt.plot(sgr_xax, melt[n,:], linewidth=1, label = iceshelves[n])
    sgr_xax = np.array([0, 1, 4, 8])
    plt.plot(sgr_xax, (melt[n,:]-melt[n,0])/np.max(np.abs(melt[n,:]-melt[n,0])), 'x', linewidth=1, linestyle = '-', label = iceshelves[n])

plt.xlabel('SGR (m^3/s)')
plt.ylabel('Integrated melf flux (GT/a)')
plt.legend(loc = 2,fontsize=8)
plt.grid()
plt.show()