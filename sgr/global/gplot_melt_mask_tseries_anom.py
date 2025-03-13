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

do_melt = 0
do_sgr = 1

tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]
sims = ["S12_control", "S12_mali", "S12_mali_x4", "S12_mali_x8"]
legsims = ["control", "sgr-x1", "sgr-x4" , "sgr-x8"]
Tmax = 50

#tseries = ["0001-0110", "0001-0110"]
#tsegment = ["clim_101-110_ts_1-110", "clim_101-110_ts_1-110"]
#sims = ["S12_control", "S12_mali"]
#legsims = ["control", "sgr-x1"]
#Tmax = 110

fyrs = f'y{Tmax}'

#iceshelves = ["Thwaites", 'Getz']
#iceshelves = ["Antarctica", "Amery", "RoiB", "Munin", "Nivl", "Fimbul", "Ekstrom", "Thwaites", "Pine_Island", "Getz", "Totten", "MoscowU", "George_VI", "Stange", "Abbot", "Larsen_C", "Filchner-Ronne", "Ross", "Riiser-Larsen/Brunt", "Antarctica-Amery"]
iceshelves = ["Antarctica", "Antarctica-Amery", "George_VI", "Stange", "Abbot", "Pine_Island", "Thwaites", "Getz", "NSS", "Ross", "Totten", "MoscowU", "Shackleton", "West", "Amery", "RoiB", "Munin", "Nivl", "Fimbul", "Ekstrom", "Riiser-Larsen/Brunt", "Filchner-Ronne","Larsen_C"]

iam, areaCell, isz = gmask_is.get_mask(iceshelves)
print(iam.shape)
iis = iam[0,:]
print(iis)

if do_sgr == 1:
    sgr_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/sgr_files_IV/MALI_20240327/DSGR.massFlux.MALI.out2055.SOwISC12to60E2r4.20240328.nc'
    ds = xarray.open_dataset(sgr_file)
    ds.load()
    sgr = np.squeeze(ds.subglacialRunoffFlux.data)

    sgr_vol_flux = np.zeros((len(iceshelves)))
    is_area = np.zeros((len(iceshelves)))
    for n in range(len(iceshelves)):
        print(iceshelves[n])
        iis = iam[n, :]
        sgr_mass_flux = np.nansum(sgr[iis] * areaCell[iis], axis=0)
        is_area[n] = np.nansum(areaCell[iis], axis=0)
        sgr_vol_flux[n] = sgr_mass_flux*secPerYear #kg per year
        print(sgr_vol_flux[n])

    print('save to netcdf')
    # Create an xarray DataArray from the NumPy array
    da = xarray.DataArray(sgr_vol_flux, dims=("iceshelves"),
                                    coords={'iceshelves': iceshelves},
                                    name='SgrMassFlux')  # Set the variable name

    #da.to_netcdf(f"/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/SgrMassFlux.nc")
    # Create an xarray DataArray from the NumPy array
    db = xarray.DataArray(is_area, dims=("iceshelves"),
                                    coords={'iceshelves': iceshelves},
                                    name='IceShelfArea')  # Set the variable name

    ds = xarray.Dataset({
        'SgrMassFlux': da,
        'IceShelfArea': db
    })

    ds.to_netcdf(f"/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/SgrMassFlux.nc")



if do_melt == 1:

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
        print('melt assigned')

        for ind_r in range(len(iceshelves)):
            print(iceshelves[ind_r])
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

    print('save to netcdf')
    # Create an xarray DataArray from the NumPy array
    da = xarray.DataArray(MeltFlux, dims=("time", "sims", "iceshelves"),
                                    coords={'time': Time,
                                        'sims': sims,
                                        'iceshelves': iceshelves},
                                    name='MeltFlux')  # Set the variable name

    da.to_netcdf(f"/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/melt_flux_tseries_{fyrs}.nc")


    print('plot')
    #fHeight = 4
    #fWidth = 7
    #plt.figure(figsize=(fWidth, fHeight))
    #plt.clf()
    clr = ["brown", "orange", "deepskyblue" , "black"]
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

