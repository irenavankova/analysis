#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os
import gmask_is
from scipy.optimize import curve_fit

do_calc_local = 0

# Get subglacial discharge per ice shelf
d2y = 365
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/global/Scaling'

iceshelves = ["George_VI", "Pine_Island", "Thwaites", "Totten", "Fimbul", "Ross", "Amery", "Filchner-Ronne", "Larsen_C", "Getz"]

nmnow = 'All'
opt_save = 0
opt_fits = 1

sgr_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/sgr_files_IV/MALI_20240327/DSGR.massFlux.MALI.out2055.SOwISC12to60E2r4.20240328.nc'
ds = xarray.open_dataset(sgr_file)
ds.load()
sgr = np.squeeze(ds.subglacialRunoffFlux.data)

p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'

dsMesh = xarray.open_dataset(p_file)
dsMesh = dsMesh[['maxLevelCell']]
dsMesh.load()
maxLevelCell = dsMesh.maxLevelCell.data - 1

iam, areaCell, isz = gmask_is.get_mask(iceshelves)

sgr_vol_flux = np.zeros((len(iceshelves)))
for n in range(len(iceshelves)):
    iis = iam[n, :]
    sgr_mass_flux = np.nansum(sgr[iis] * areaCell[iis], axis=0)
    sgr_vol_flux[n] = sgr_mass_flux/rho_fw
    print(sgr_vol_flux[n])

# Get melt rates
sims = ["S12_control", "S12_mali", "S12_mali_x4" , "S12_mali_x8"]
#sims = ["S12_control", "S12_mali", "S12_mali_x4" , "S12_mali_x10"]


tseg = 9
#c1 = np.array([0, 0, 0, 0])+ts1
c1 = np.array([0, 0, 0, 0])+41
c2 = c1 + tseg
t2 = c1 + tseg

tseries = []
tsegment = []
tclim = []

for n in range(len(c1)):
    tseries.append(f"0001-{t2[n]:04}")
    tsegment.append(f"clim_{c1[n]}-{c2[n]}_ts_1-{t2[n]}")
    tclim.append(f"{c1[n]:04}01_{c2[n]:04}12")

#tseries = ["0001-0110", "0001-0110", "0001-0050", "0001-0050"]
#tsegment = [f'clim_{c1[0]}-110_ts_1-110', "clim_101-110_ts_1-110", "clim_41-50_ts_1-50", "clim_41-50_ts_1-50"]

melt = np.zeros((len(iceshelves), len(sims)))
mr_vol_flux = np.zeros((len(iceshelves), len(sims)))
tmax = np.zeros((len(iceshelves), len(sims)))

#for n in range(1):
for n in range(len(iceshelves)):
    shelf_name = iceshelves[n]
    if shelf_name == "PineIsland":
        shelf_name = "Pine_Island"

    #for s in range(1):
    for s in range(len(sims)):
        #print(shelf_name)
        iis = iam[n, :]

        c_file = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/mpaso_ANN_{tclim[s]}_climo.nc'
        dsOut = xarray.open_dataset(c_file)
        dsOut = dsOut[['timeMonthly_avg_landIceFreshwaterFluxTotal','timeMonthly_avg_activeTracers_temperature']]
        dsOut.load()
        mr_clim = np.squeeze(dsOut.timeMonthly_avg_landIceFreshwaterFluxTotal.data)
        temp_clim = np.squeeze(dsOut.timeMonthly_avg_activeTracers_temperature.data)
        ix = np.where(iis)[0]
        '''
        print(iis.shape)
        print(ix.shape)
        print(temp_clim[ix[0],:])
        print(temp_clim[ix[0], maxLevelCell[ix[0]]])
        print(temp_clim[ix[2], :])
        print(temp_clim[ix[2], maxLevelCell[ix[2]]])
        print(temp_clim)
        print('BREAK')
        
        fHeight = 5
        fWidth = 8
        plt.figure(figsize=(fWidth, fHeight))
        plt.hist(temp_clim)
        plt.show()
        '''
        # Get max sea floor temperature
        temp_clim = temp_clim[ix,maxLevelCell[ix]]
        #tmax[n, s] = np.max(temp_clim, axis=0)
        tmax[n, s] = np.mean(temp_clim, axis=0)
        # Get melt rate
        mr_mass_flux = np.nansum(mr_clim[iis] * areaCell[iis], axis=0)
        #mr_vol_flux[n,s] = mr_mass_flux / rho_fw * secPerYear
        mr_vol_flux[n, s] = mr_mass_flux * secPerYear / 10.0**12


def f_n_pow_1_3(x, a):
    return a * np.power(x, 1/3)

def f_n_pow_2_3(x, a):
    return a * np.power(x, 2/3)

xfit = np.linspace(0,1200,100)

fHeight = 5
fWidth = 8
clr = 'rkgbmcyrkgb'
#smb = 'ooooo^^^sss'
smb = 'ppppp^^^sss'
sgr_mat = np.zeros((len(iceshelves), len(sims)))
for n in range(len(iceshelves)):
    #sgr_xax = np.array([0, 1, 4, 8]) * sgr_vol_flux[n]
    #plt.plot(sgr_xax, melt[n,:], linewidth=1, label = iceshelves[n])
    sgr_xax = np.array([0, 1, 4, 8]) * sgr_vol_flux[n]
    sgr_mat[n,:] = sgr_xax

    '''
    if opt_fits == 1:
        plt.plot(sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0], f'{clr[n]}{smb[n]}', linewidth=1, linestyle=':')
        popt, pcov = curve_fit(f_n_pow_1_3, sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0], p0=10 ** 10)
        plt.plot(xfit, f_n_pow_1_3(xfit, *popt), f'{clr[n]}--', linewidth=1)
        popt, pcov = curve_fit(f_n_pow_2_3, sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0])
        plt.plot(xfit, f_n_pow_2_3(xfit, *popt), f'{clr[n]}-', linewidth=1)
    else:
        plt.plot(sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0], f'{clr[n]}{smb[n]}', linewidth=1, linestyle=':')
    
    fig, (ax1, ax2) = plt.subplots(2, 1)  # 2 rows, 1 column
    ax1.plot(sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0], f'{clr[n]}{smb[n]}', linewidth=1, linestyle=':')
    ax2.plot(sgr_xax, tmax[n, :], f'{clr[n]}{smb[n]}', linewidth=1, linestyle=':')
    plt.show()
    '''

#Plot of T vs SGR points
fig1, ax1 = plt.subplots()
#svec = sgr_mat.flatten()
#tvec = tmax.flatten()
#ax1.scatter(svec, tvec)
for n in range(len(iceshelves)):
    ax1.scatter(sgr_mat[n,:], tmax[n,:],label = iceshelves[n])

plt.legend(loc = 'upper right' ,fontsize=8)
plt.show()

'''
plt.xlabel('SGR (m^3/s)')
plt.ylabel('Integrated melf flux (GT/a)')
plt.legend(loc = 1,fontsize=8)
plt.grid()
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/{nmnow}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()
'''

#plt.figure(figsize=(fWidth, fHeight))
#for n in range(len(iceshelves)):
    #plt.plot(sgr_xax, mr_vol_flux[n, :] - mr_vol_flux[n, 0], f'{clr[n]}{smb[n]}', linewidth=1, linestyle=':')


'''
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
if opt_save != 1:
    plt.show()
'''
