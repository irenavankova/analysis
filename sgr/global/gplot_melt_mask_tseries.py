#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
from scipy import signal
import gmask_is

opt_save = 1
dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/global/Tseries'
Lwide = 0.75

d2y = 365
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

ny = 110
tsegment = 'clim_101-110_ts_1-110'
tseries = f'0001-0110'
#sims = ["S12_mali"]
sims = ["S12_control", "S12_mali"]
#ny = 32
#tsegment = 'clim_23-26_ts_1-32'
#tseries = f'0001-0032'
#sims = ["S12_control"]

Time = (np.arange(12*ny)+1)/12

#is_list = ["Amery" , "RoiB", "Munin", "Nivl", "Fimbul", "Ekstrom"]
#is_list = ["Amery" , "RoiB", "Munin", "Fimbul"]

casenum = 'Amundsen'

if casenum == 'Amundsen':
    is_list = ["Pine Island" , "Thwaites", "Getz"]
elif casenum == 'DroningMaud':
    is_list = ["Amery", "RoiB", "Munin", "Fimbul"]


iam, areaCell = gmask_is.get_mask(is_list)
print(iam.shape)
iis = iam[0,:]
print(iis)

for s in range(len(sims)):
    #Load melt
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment}/mpasTimeSeriesOcean.nc'
    dsOut = xarray.open_dataset(p_file)
    dsOut = dsOut[['timeMonthly_avg_landIceFreshwaterFluxTotal']]
    dsOut.load()
    melt = np.squeeze(dsOut.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
    melt = melt / rho_fw * secPerYear
    if s == 0:
        melt_flux = np.zeros((len(Time), len(is_list), len(sims)))

    # Calculate melt flux timeseries for each ice shelf of interest
    for n in range(len(is_list)):
        iis = iam[n,:]
        print(np.nansum(melt[:,iis], axis=1))
        melt_flux[:,n,s] = np.nansum(melt[:,iis] * areaCell[iis], axis=1)

fHeight = 5
fWidth = 8

plt.figure(figsize=(fWidth, fHeight))

for s in range(len(sims)):
    if s == 0:
        smb = '-'
    else:
        smb = ':'

    for n in range(len(is_list)):
        plt.plot(Time, np.squeeze(melt_flux[:,n,s]),smb, linewidth=Lwide, label = is_list[n])

fsize = 8
plt.legend(loc = 2,fontsize=fsize)
plt.xlim([19, Time[-1]])
plt.ylabel('Melt-flux (m$^3$/a)',fontsize=fsize)
plt.xlabel('Time (a)',fontsize=fsize)
plt.rcParams.update({'font.size': fsize})
plt.tick_params(axis='both', labelsize=fsize)
plt.grid()

if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/MFlux_{casenum}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

plt.figure(figsize=(fWidth, fHeight))

smb = '-'
for n in range(len(is_list)):
    plt.plot(Time, np.squeeze(melt_flux[:,n,1])-np.squeeze(melt_flux[:,n,0]),smb, linewidth=Lwide, label = is_list[n])

fsize = 8
plt.legend(loc = 2,fontsize=fsize)
plt.xlim([19, Time[-1]])
plt.ylabel('Melt-flux anomaly (m$^3$/a)',fontsize=fsize)
plt.xlabel('Time (a)',fontsize=fsize)
plt.rcParams.update({'font.size': fsize})
plt.tick_params(axis='both', labelsize=fsize)
plt.grid()
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/MFlux_anom_{casenum}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

if casenum == 'Amundsen':
    ind_T = (Time > 20)
    for k in range(2):
        n = 0
        xp = np.squeeze(melt_flux[ind_T,n,1])-np.squeeze(melt_flux[ind_T,n,0])
        xp = signal.detrend(xp)
        if k == 0:
            n = 2
            ISlagged = '-Getz'
            mp = -1
        else:
            n = 1
            ISlagged = 'Thwaites'
            mp = 1
        xg = mp*(np.squeeze(melt_flux[ind_T,n,1])-np.squeeze(melt_flux[ind_T,n,0]))
        xg = signal.detrend(xg)

        # Calculate correlation
        correlation = np.correlate(xp, xg, mode='full')
        correlation = correlation/(np.sqrt(np.sum(np.abs(xg) ** 2) * np.sum(np.abs(xp) ** 2)))

        # Get lag values
        lags = signal.correlation_lags(xp.size, xg.size, mode='full')

        # Find lag with maximum correlation
        max_correlation_index = np.argmax(correlation)
        max_lag = lags[max_correlation_index]

        print("Correlation:", correlation)
        print("Lags:", lags)
        print("Maximum correlation:", correlation[max_correlation_index])
        print("Lag at maximum correlation:", max_lag)

        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(Time[ind_T],xp, linewidth=Lwide, label = 'Pine Island')
        plt.plot(Time[ind_T],xg,':', linewidth=Lwide, label = f'{ISlagged}')
        plt.plot(Time[ind_T]+max_lag/12,xg, linewidth=Lwide, label = f'{ISlagged} lagged')
        plt.legend(loc = 2,fontsize=fsize)
        plt.ylabel('Detrended melt-flux anomaly (m$^3$/a)',fontsize=fsize)
        plt.xlabel('Time (a)',fontsize=fsize)
        plt.rcParams.update({'font.size': fsize})
        plt.tick_params(axis='both', labelsize=fsize)
        #plt.text(.05, .95, f'lag = {max_lag}', ha='right', va='top', transform=ax.transAxes)
        plt.grid()
        plt.title(f'lag = {max_lag} months, xcor = {np.round(correlation[max_correlation_index],decimals=2)}', fontsize=fsize)

        if opt_save == 1:
            plt.savefig(f'{dir_fig_save}/MFlux_anom_PIG_{ISlagged}_lag.png', bbox_inches='tight', dpi=600)
        else:
            plt.show()

if casenum == 'Amundsen':
    ind_T = (Time > 20)
    for k in range(2):
        n = 0
        xp = np.squeeze(melt_flux[ind_T,n,0])
        xp = signal.detrend(xp)
        if k == 0:
            n = 2
            ISlagged = 'Getz'
            mp = 1
        else:
            n = 1
            ISlagged = 'Thwaites'
            mp = 1
        xg = mp*np.squeeze(melt_flux[ind_T,n,0])
        xg = signal.detrend(xg)

        # Calculate correlation
        correlation = np.correlate(xp, xg, mode='full')
        correlation = correlation/(np.sqrt(np.sum(np.abs(xg) ** 2) * np.sum(np.abs(xp) ** 2)))

        # Get lag values
        lags = signal.correlation_lags(xp.size, xg.size, mode='full')

        # Find lag with maximum correlation
        max_correlation_index = np.argmax(correlation)
        max_lag = lags[max_correlation_index]

        print("Correlation:", correlation)
        print("Maximum correlation:", correlation[max_correlation_index])
        print("Lag at maximum correlation:", max_lag)

        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(Time[ind_T],xp, linewidth=Lwide, label = 'Pine Island')
        plt.plot(Time[ind_T],xg,':', linewidth=Lwide, label = f'{ISlagged}')
        plt.plot(Time[ind_T]+max_lag/12,xg, linewidth=Lwide, label = f'{ISlagged} lagged')
        plt.legend(loc = 2,fontsize=fsize)
        plt.ylabel('Detrended melt-flux control (m$^3$/a)',fontsize=fsize)
        plt.xlabel('Time (a)',fontsize=fsize)
        plt.rcParams.update({'font.size': fsize})
        plt.tick_params(axis='both', labelsize=fsize)
        #plt.text(.05, .95, f'lag = {max_lag}', ha='right', va='top', transform=ax.transAxes)
        plt.grid()
        plt.title(f'lag = {max_lag} months, xcor = {np.round(correlation[max_correlation_index],decimals=2)}', fontsize=fsize)

        if opt_save == 1:
            plt.savefig(f'{dir_fig_save}/MFlux_control_PIG_{ISlagged}_lag.png', bbox_inches='tight', dpi=600)
        else:
            plt.show()




