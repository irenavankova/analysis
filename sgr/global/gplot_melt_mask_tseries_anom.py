#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import gmask_is
import cftime


d2y = 365
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

#tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
#tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
#sims = ["S12_control", "S12_mali", "S12_mali_x4", "S12_mali_x8"]
#legsims = ["control", "sgr-x1", "sgr-x4" , "sgr-x8"]

tseries = ["0001-0110", "0001-0110"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110"]
sims = ["S12_control", "S12_mali"]
legsims = ["control", "sgr-x1"]

iceshelves = ["Thwaites"]
iam, areaCell, isz = gmask_is.get_mask(iceshelves)
print(iam.shape)
iis = iam[0,:]
print(iis)

fHeight = 4
fWidth = 7
plt.figure(figsize=(fWidth, fHeight))
plt.clf()
clr = ["brown", "orange", "deepskyblue" , "black"]

Tmax = 50
Tmin = 20

s = 1
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/LIFW.nc'
dsOut = xarray.open_dataset(p_file)
Time = np.arange(1,Tmax*12+1,1)/12
print(len(Time))

MeltFlux = np.zeros((len(Time), len(sims), len(iceshelves)))

for s in range(len(sims)):
    print(s)
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/LIFW.nc'
    dsOut = xarray.open_dataset(p_file)
    dsOut = dsOut[['timeMonthly_avg_landIceFreshwaterFluxTotal']]
    dsOut.load()
    print('melt loaded')
    melt = np.squeeze(dsOut.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
    melt = melt / kgingt * secPerYear
    print('melt assighen')

    for ind_r in range(len(iceshelves)):
        print(ind_r)
        iis = iam[ind_r,:]
        print(len(iis))
        ix = np.where(iis)[0]
        print(len(ix))
        print(ix)
        #dsIs = dsOut.sel(Time=dsOut.Time[dsOut.Time <= d2y * Tmax], nCell=[iis])
        #dsIs = dsOut.sel(nCells=[iis])
        #dsIs = dsOut.isel(nCells=ix)
        #dsIs.load()
        #print(dsIs)
        #print('MeltFluxNow')
        #MeltFluxNow = np.squeeze(dsIs.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
        print('MeltFluxNow loaded')
        MeltFluxNow = melt[0:len(Time),ix]
        #MeltFluxNow = MeltFluxNow / rho_fw * secPerYear

        print(MeltFluxNow.shape)
        print(areaCell[ix].shape)
        print('before sum')
        #MeltFlux[:, s, ind_r] = np.nansum(MeltFluxNow, axis=1)
        MeltFlux[:, s, ind_r] = np.nansum(MeltFluxNow * areaCell[ix], axis=1)
        print('after sum')

print('rho')
#MeltFlux = MeltFlux / rho_fw * secPerYear

print('plot')
for ind_r in range(len(iceshelves)):

    for s in range(len(sims)):

        if s == 0:
            #smb = 'k:'
            smb = '-'
            ref = np.squeeze(MeltFlux[:, s, ind_r])
        else:
            smb = '-'

        print(np.squeeze(MeltFlux[:, s, ind_r]))
        plt.plot(Time, np.squeeze(MeltFlux[:, s, ind_r]) ,smb, color = clr[s] ,linewidth=1.5, label = legsims[s])
        #ind_t = np.where((Time >= 25) & (Time <= 50))
        #mmean = np.mean(MeltFluxNow[ind_t, ind_r])
        #print(mmean)

    plt.xlim([Tmin, Tmax])
    plt.xlabel('Time (a)')
    plt.ylabel('Integrated melf flux (GT/a)')
    plt.legend(loc = 2,fontsize=12)
    #plt.title(regionNames[ind_r])
    plt.show()

# Create an xarray DataArray from the NumPy array
da = xarray.DataArray(MeltFlux, dims=("time", "sims", "iceshelves"),
                                coords={'time': Time,
                                    'sims': sims,
                                    'iceshelves': iceshelves},
                                name='MeltFlux')  # Set the variable name

da.to_netcdf("/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/my_array.nc")