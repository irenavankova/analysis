#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

opt_save = 0
d2y = 365


tsegment = 'clim_23-26_ts_1-32'
tseries = '0001-0032'
sims = ["S12_control", "S12_mali" , "S12_mali_x4" , "S12_mali_x10"]
#sims = ["S12_mali_x4"]

'''
tsegment = 'clim_101-110_ts_1-110'
tseries = f'0001-0110'
sims = ["S12_control", "S12_mali"]
'''

for s in range(len(sims)):
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment}/iceShelfFluxes_{tseries}.nc'

    dsOut = xarray.open_dataset(p_file)
    dsOut.load()
    MeltFluxNow = np.squeeze(dsOut.integratedMeltFlux.data)

    if s == 0:
        regionNames = dsOut.regionNames.data
        regionNames = np.squeeze(regionNames[0,:])
        ind_x = np.where((regionNames == "Totten"))
        print(regionNames[ind_x])
        Time = np.squeeze(dsOut.Time.data)/d2y
        MeltFlux = np.zeros((len(Time), len(sims)))

    MeltFlux[:, s] = np.squeeze(MeltFluxNow[:,ind_x])

fHeight = 5
fWidth = fHeight


plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)):
    plt.plot(Time, np.squeeze(MeltFlux[:,s]))

plt.xlim([19, Time[-1]])
plt.show()
'''
plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)-1):
    plt.plot(Time, np.squeeze(MeltFlux[:,s+1])-np.squeeze(MeltFlux[:,0]))

plt.xlim([19, Time[-1]])
plt.show()
'''
ind_t = np.where((Time > 22))

plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)-1):
    plt.plot(np.squeeze(MeltFlux[ind_t,0]), np.squeeze(MeltFlux[ind_t,s+1])-np.squeeze(MeltFlux[ind_t,0]),'.', label = sims[s+1])

plt.legend(loc = 2)
#plt.xlim([19, Time[-1]])
plt.show()



