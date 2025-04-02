#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import cmocean

opt_save = 1

rho_fw = 1000.
kgInGt = 1e12

Nmali = [0, 1, 4, 8]
clr = ["brown", "orange", "deepskyblue" , "black"]

cdir = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/post_derived/sgr_fw/fw_tseries_M'

ncols = 2
nrows = 1

fHeight = 6
fWidth = 16
cm = 1/2.54

fig, ax = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.1*cm,wspace=0.1*cm)  # Increase vertical spacing
fsize = 8

for k in range(2):
    for j in range(len(Nmali)):
        p_file = f'{cdir}{Nmali[j]}.nc'
        ds = xr.open_dataset(p_file)
        ds.load()

        if k == 0:
            ttle = 'Land ice melt'
            icav = np.squeeze(ds['VFWLIicav'].values)
            ishelf = np.squeeze(ds['VFWLIishelf'].values)
            is60 = np.squeeze(ds['VFWLI1000is60'].values)
        else:
            ttle = 'Subglacial discharge'
            icav = np.squeeze(ds['VFWSRicav'].values)
            ishelf = np.squeeze(ds['VFWSRishelf'].values)
            is60 = np.squeeze(ds['VFWSR1000is60'].values)

        # m^3 to Gt
        icav = icav * rho_fw / kgInGt
        ishelf = ishelf * rho_fw / kgInGt
        is60 = is60 * rho_fw / kgInGt

        if (Nmali[j] == 1) & (k == 0):
            print(Nmali[j])
            print(k)
            ymax = np.nanmax(is60)

        time = (np.squeeze(ds['time'].values)+1)/12

        ax[k].plot(time, icav, color = clr[j] ,linewidth=0.75, linestyle='-')#, label=f'icav{Nmali[j]}')
        ax[k].plot(time, ishelf, color=clr[j], linewidth=0.75, linestyle='--')#, label=f'ishelf{Nmali[j]}')
        ax[k].plot(time, is60, color=clr[j], linewidth=0.75, linestyle=':')#, label=f'is60{Nmali[j]}')

    #ax[k].legend(loc=2, prop={'size': fsize})
    ax[k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    ax[k].set_title(f'{ttle}', fontsize=fsize)
    if k == 0:
        #ax[k].set_ylabel(f'Freshwater volume (m$^3$)', fontsize=fsize)
        ax[k].set_ylabel(f'Freshwater mass (GT)', fontsize=fsize)
    else:
        plt.setp(ax[k].get_yticklabels(), visible=False)

    ax[k].set_xlabel('Time (a)', fontsize=fsize)
    ax[k].tick_params(axis='both', labelsize=fsize)
    # ax[k].rcParams.update({'font.size': 8})
    ax[k].set_xlim([0, 90])
    ax[k].set_ylim([0, ymax])




if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/melt_tseries/FW_tseries.png', bbox_inches='tight', dpi=300)
else:
    plt.show()