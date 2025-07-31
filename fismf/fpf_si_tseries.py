#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt
import matplotlib.ticker as ticker


opt_save = 1
fname = 'pANS'
ff = 1
var_plot = 'siv'
if var_plot == 'sic':
    fac = 1e6
else:
    fac = 1e9
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

#opt_asc = 'asc_transp'
opt_asc = ''

fname_hist = f'si_reg_hist_ave_tseries.nc'
fname_pismf = f'si_reg_pismf_ave_tseries.nc'
fname_fismf = f'si_reg_fismf_ave_tseries.nc'

ds_hist = xr.open_dataset(f'{fpath}{fname_hist}')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}')

# Extract the DataArray (assuming variable is 'lifw')
Fh = ds_hist[var_plot]  # dims: (nTransects, Time)
n_time = len(ds_hist['Time'])
Fht = np.arange(0,n_time)/12 + 1951

Fp = ds_pismf[var_plot]  # dims: (nTransects, Time)
Ff = ds_fismf[var_plot]  # dims: (nTransects, Time)
n_time = len(ds_pismf['Time'])
Fpt = np.arange(0,n_time)/12 + 2015

fs = 1/(Fht[1]-Fht[0])
fc = 1/((Fht[1]-Fht[0])*36)
clr = ["black", "darkorange","lightskyblue", "brown", "royalblue"]
smb = ["-", "--", "--", "-", "-"]

si_max_mo = 9
si_min_mo = 2

reg_names = ["Southern Ocean", "Bellings./Amundsen", "Eastern Ross", "Western Ross", "East Antarctica", "Amery", "Dronning Maud Land", "Weddell"]
abc = 'abcdefghijklmnop'
ncols = 2
nrows = int(np.ceil((len(ds_pismf['region']))) / ncols)

fHeight = 16
fWidth = 16
cm = 1/2.54
Lwide = 0.5
Lwide2 = 1.0

fig, axes = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.3)  # Increase vertical spacing
plt.subplots_adjust(wspace=0.2)  # Increase vertical spacing

j = int(0); k = int(0)

for r, region in enumerate(Fh.region.values):

    Fh_reg = Fh.sel(region=region) / fac
    Ff_reg = Ff.sel(region=region) / fac
    Fp_reg = Fp.sel(region=region) / fac

    axes[j, k].plot(Fht, Fh_reg, '--', color="yellowgreen", linewidth=Lwide)
    axes[j, k].plot(Fpt, Fp_reg, '--', color="darkorange", linewidth=Lwide)
    axes[j, k].plot(Fpt, Ff_reg, '--', color="lightskyblue", linewidth=Lwide)


    if ff == 1:
        Fh_filt = Fh_reg[si_max_mo::12]
        Fp_filt = Fp_reg[si_max_mo::12]
        Ff_filt = Ff_reg[si_max_mo::12]
        Fht_filt = Fht[si_max_mo::12]
        Fpt_filt = Fpt[si_max_mo::12]

        axes[j, k].plot(Fht_filt, Fh_filt, '-', color="darkolivegreen", linewidth=Lwide2, label='HIST')
        axes[j, k].plot(Fpt_filt, Fp_filt, '-', color="brown", linewidth=Lwide2, label='EVOMELT')
        axes[j, k].plot(Fpt_filt, Ff_filt, '-', color="royalblue", linewidth=Lwide2, label='FIXMELT')

        Fh_filt = Fh_reg[si_min_mo::12]
        Fp_filt = Fp_reg[si_min_mo::12]
        Ff_filt = Ff_reg[si_min_mo::12]
        Fht_filt = Fht[si_min_mo::12]
        Fpt_filt = Fpt[si_min_mo::12]

        axes[j, k].plot(Fht_filt, Fh_filt, '-', color="darkolivegreen", linewidth=Lwide2)
        axes[j, k].plot(Fpt_filt, Fp_filt, '-', color="brown", linewidth=Lwide2)
        axes[j, k].plot(Fpt_filt, Ff_filt, '-', color="royalblue", linewidth=Lwide2)

    axes[j, k].grid(True)

    fsize = 8
    axes[j, k].set_xlim(np.array([np.min(Fht), np.max(Fpt)]))
    axes[j, k].set_title(f'{abc[r]}) {reg_names[r]}', fontsize=fsize - 1)
    # axes[j, k].autoscale(enable=True, axis='both', tight=True)
    axes[j, k].tick_params(axis='both', labelsize=fsize)

    axes[j, k].grid(which='major', linestyle=':', linewidth='0.5', color='gray')
    axes[j, k].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))

    if k == 0:
        if var_plot == 'sic':
            axes[j, k].set_ylabel('Sea Ice Area (km^2)', fontsize=fsize)
        else:
            axes[j, k].set_ylabel('SI Vol (km$^3$)', fontsize=fsize)

    if j == nrows - 1:
        axes[j, k].set_xlabel('Time (a)', fontsize=fsize)
    else:
        plt.setp(axes[j, k].get_xticklabels(), visible=False)

    if r == 0:
        axes[j, k].legend()
        axes[j, k].legend(fontsize=6, loc='lower left')

    #if r == 0:
    #    fig.delaxes(axes[0, 1])
    #    k = 0
    #    j = 1
    #else:
    k = k + 1
    if k > ncols - 1:
        k = 0
        j = j + 1

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/SeaIceTseries/SI_tseries_{fname}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()