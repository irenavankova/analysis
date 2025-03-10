#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

d2y = 365
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
#tseries = ["0001-0050", "0001-0050", "0001-0050", "0001-0050"]
#tsegment = ["clim_41-50_ts_1-50", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
sims = ["S12_control", "S12_mali", "S12_mali_x4", "S12_mali_x8"]
legsims = ["control", "sgr-x1", "sgr-x4" , "sgr-x8"]

fHeight = 4
fWidth = 7
plt.figure(figsize=(fWidth, fHeight))
plt.clf()
clr = ["brown", "orange", "deepskyblue" , "black"]

s = 3
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/iceShelfFluxes_{tseries[s]}.nc'
dsOut = xarray.open_dataset(p_file)
dsOut.load()
regionNames = dsOut.regionNames.data
regionNames = np.squeeze(regionNames[0,:])
Time = np.squeeze(dsOut.Time.data) / d2y

MeltFlux = np.zeros((len(Time), len(sims), len(regionNames)))

for s in range(len(sims)):
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/iceShelfFluxes_{tseries[s]}.nc'
    dsOut = xarray.open_dataset(p_file)
    #dsOut.load()
    dsOut = dsOut.sel(Time=dsOut.Time[dsOut.Time <= d2y * 50])
    MeltFluxNow = np.squeeze(dsOut.integratedMeltFlux.data)

    for ind_r in range(len(regionNames)):
        MeltFlux[:, s, ind_r] = np.squeeze(MeltFluxNow[:, ind_r])
        #MeltFlux[:, s, ind_r] = np.squeeze(MeltFluxNow[0:len(Time),ind_r])

for ind_r in range(len(regionNames)):

    for s in range(len(sims)):

        if s == 0:
            #smb = 'k:'
            smb = '-'
            ref = np.squeeze(MeltFlux[:, s, ind_r])
        else:
            smb = '-'

        plt.plot(Time, np.squeeze(MeltFlux[:, s, ind_r])/ref ,smb, color = clr[s] ,linewidth=1.5, label = legsims[s])
        #ind_t = np.where((Time >= 25) & (Time <= 50))
        #mmean = np.mean(MeltFluxNow[ind_t, ind_r])
        #print(mmean)

    plt.xlim([20, 50])
    plt.xlabel('Time (a)')
    plt.ylabel('Integrated melf flux (GT/a)')
    plt.legend(loc = 2,fontsize=12)
    plt.title(regionNames[ind_r])
    plt.show()




