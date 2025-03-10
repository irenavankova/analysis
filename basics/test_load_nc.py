#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os


p_file = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/my_array.nc'
dsOut = xarray.open_dataset(p_file)

MeltFlux = np.array(dsOut.MeltFlux.data)
Time = np.array(dsOut.time.data)
iceshelves = np.array(dsOut.iceshelves.data)
sims = np.array(dsOut.sims.data)

print(iceshelves)
print(iceshelves[0])

Tmax = 50
Tmin = 20


fHeight = 4
fWidth = 7
plt.figure(figsize=(fWidth, fHeight))
plt.clf()
clr = ["brown", "orange", "deepskyblue" , "black"]

#np.shape(H)[1]
#for ind_r in range(np.shape(MeltFlux)[2]):
for ind_r in range(1):

    for s in range(len(sims)):

        if s == 0:
            #smb = 'k:'
            smb = '-'
            ref = np.squeeze(MeltFlux[:, s, ind_r])
        else:
            smb = '-'

        print(np.squeeze(MeltFlux[:, s, ind_r]))
        plt.plot(Time, np.squeeze(MeltFlux[:, s, ind_r]) ,smb, color = clr[s] ,linewidth=1.5, label = sims[s])
        #ind_t = np.where((Time >= 25) & (Time <= 50))
        #mmean = np.mean(MeltFluxNow[ind_t, ind_r])
        #print(mmean)

    plt.xlim([Tmin, Tmax])
    plt.xlabel('Time (a)')
    plt.ylabel('Integrated melf flux (GT/a)')
    plt.legend(loc = 2,fontsize=12)
    plt.title(iceshelves[ind_r])
    plt.show()