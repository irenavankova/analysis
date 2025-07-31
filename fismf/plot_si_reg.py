#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt

opt_save = 0
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

for region in Fh.region.values:
    plt.figure(figsize=(6, 3))

    Fh_reg = Fh.sel(region=region) / fac
    Ff_reg = Ff.sel(region=region) / fac
    Fp_reg = Fp.sel(region=region) / fac

    Lwide = 1.0

    plt.plot(Fht, Fh_reg, '--', color="yellowgreen", linewidth=Lwide)
    plt.plot(Fpt, Fp_reg, '--', color="darkorange", linewidth=Lwide)
    plt.plot(Fpt, Ff_reg, '--', color="lightskyblue", linewidth=Lwide)


    if ff == 1:
        Fh_filt = Fh_reg[si_max_mo::12]
        Fp_filt = Fp_reg[si_max_mo::12]
        Ff_filt = Ff_reg[si_max_mo::12]
        Fht_filt = Fht[si_max_mo::12]
        Fpt_filt = Fpt[si_max_mo::12]

        plt.plot(Fht_filt, Fh_filt, '-', color="darkolivegreen", linewidth=Lwide*2, label='HIST')
        plt.plot(Fpt_filt, Fp_filt, '-', color="brown", linewidth=Lwide*2, label='EVOMELT')
        plt.plot(Fpt_filt, Ff_filt, '-', color="royalblue", linewidth=Lwide*2, label='FIXMELT')

        Fh_filt = Fh_reg[si_min_mo::12]
        Fp_filt = Fp_reg[si_min_mo::12]
        Ff_filt = Ff_reg[si_min_mo::12]
        Fht_filt = Fht[si_min_mo::12]
        Fpt_filt = Fpt[si_min_mo::12]

        plt.plot(Fht_filt, Fh_filt, '-', color="darkolivegreen", linewidth=Lwide*2)
        plt.plot(Fpt_filt, Fp_filt, '-', color="brown", linewidth=Lwide*2)
        plt.plot(Fpt_filt, Ff_filt, '-', color="royalblue", linewidth=Lwide*2)

    fsize = 8
    plt.xlim(np.array([np.min(Fht), np.max(Fpt)]))
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tick_params(axis='both', labelsize=fsize)

    plt.title(region,fontsize=fsize)
    plt.xlabel('Time (a)', fontsize=fsize)
    if var_plot == 'sic':
        plt.ylabel('Sea Ice Area (km^2)', fontsize=fsize)
    else:
        plt.ylabel('Sea Ice Volume (km^3)', fontsize=fsize)

    plt.grid(True)
    plt.legend()
    plt.legend(fontsize=6)

    if opt_save == 1:
        plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/Transport/Transp_{regname.replace(" ", "_")}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()
