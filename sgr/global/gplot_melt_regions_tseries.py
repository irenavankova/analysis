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

shelf_name = "Pine_Island"
#tsegment = 'clim_23-26_ts_1-32'
#tseries = '0001-0032'
#sims = ["S12_control", "S12_mali" , "S12_mali_x4" , "S12_mali_x10"]
#sims = ["S12_mali_x4"]

tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
sims = ["S12_control", "S12_mali", "S12_mali_x4" , "S12_mali_x8"]

fHeight = 5
fWidth = 10
plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)):
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/iceShelfFluxes_{tseries[s]}.nc'

    dsOut = xarray.open_dataset(p_file)
    dsOut.load()
    MeltFluxNow = np.squeeze(dsOut.integratedMeltFlux.data)
    Time = np.squeeze(dsOut.Time.data) / d2y
    if s == 0:
        regionNames = dsOut.regionNames.data
        regionNames = np.squeeze(regionNames[0,:])
        ind_x = np.where((regionNames == shelf_name))
        print(regionNames[ind_x])
        #Time = np.squeeze(dsOut.Time.data)/d2y
        #MeltFlux = np.zeros((len(Time), len(sims)))

    #MeltFlux[:, s] = np.squeeze(MeltFluxNow[:,ind_x])
    plt.plot(Time, np.squeeze(MeltFluxNow[:, ind_x]), linewidth=1, label = sims[s])


#for s in range(len(sims)):
    #plt.plot(Time, np.squeeze(MeltFlux[:,s]))

#plt.xlim([19, Time[-1]])
plt.xlim([0, 100])
plt.title(shelf_name)
plt.xlabel('Time (a)')
plt.ylabel('Integrated melf flux (GT/a)')
plt.legend(loc = 2,fontsize=8)
plt.grid()

opt_save = 1
if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Tseries/iceshelf/{shelf_name}.png', bbox_inches='tight',
                dpi=600)
else:
    plt.show()

'''
plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)-1):
    plt.plot(Time, np.squeeze(MeltFlux[:,s+1])-np.squeeze(MeltFlux[:,0]))

plt.xlim([19, Time[-1]])
plt.show()
'''

'''
ind_t = np.where((Time > 22))

plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)-1):
    plt.plot(np.squeeze(MeltFlux[ind_t,0]), np.squeeze(MeltFlux[ind_t,s+1])-np.squeeze(MeltFlux[ind_t,0]),'.', label = sims[s+1])

plt.legend(loc = 2)
#plt.xlim([19, Time[-1]])
plt.show()
'''



