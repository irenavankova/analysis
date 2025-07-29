#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt

opt_save = 0
ff = 1
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/'

opt_asc = 'asc_transp'
#opt_asc = ''

fname_hist = f'hist/{opt_asc}/TransportTransects_1951-2014.nc'
fname_pismf = f'pismf/{opt_asc}/TransportTransects_2015-2100.nc'
#fname_fismf_701 = 'fismf_701/asc_transp/TransportTransects_2015-2100.nc'
fname_fismf = f'fismf/{opt_asc}/TransportTransects_2015-2100.nc'


ds_hist = xr.open_dataset(f'{fpath}{fname_hist}')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}')

# Extract the DataArray (assuming variable is 'lifw')
Fh = ds_hist['transport']  # dims: (nTransects, Time)
n_time = len(ds_hist['Time'])
Fht = np.arange(0,n_time)/12 + 1951
Fnames = ds_hist['transectNames'].values  # dims: (nTransects, Time)

Fp = ds_pismf['transport']  # dims: (nTransects, Time)
Ff = ds_fismf['transport']  # dims: (nTransects, Time)
n_time = len(ds_pismf['Time'])
Fpt = np.arange(0,n_time)/12 + 2015

fs = 1/(Fht[1]-Fht[0])
fc = 1/((Fht[1]-Fht[0])*36)
clr = ["black", "darkorange","lightskyblue", "brown", "royalblue"]
smb = ["-", "--", "--", "-", "-"]

for transect in Fh.nTransects.values:
    plt.figure(figsize=(6, 3))

    Fh_reg = Fh.sel(nTransects=transect)
    Ff_reg = Ff.sel(nTransects=transect)
    Fp_reg = Fp.sel(nTransects=transect)

    Lwide = 1.0

    plt.plot(Fht, Fh_reg, '--', color="yellowgreen", linewidth=Lwide)
    plt.plot(Fpt, Fp_reg, '--', color="darkorange", linewidth=Lwide)
    plt.plot(Fpt, Ff_reg, '--', color="lightskyblue", linewidth=Lwide)

    if ff == 1:
        p_reg = np.append(Fh_reg, Fp_reg, axis=0)
        f_reg = np.append(Fh_reg, Ff_reg, axis=0)

        p_reg = iv_filt.butter_filter(p_reg, fc, fs, 'low')
        f_reg = iv_filt.butter_filter(f_reg, fc, fs, 'low')

        Fh_filt = p_reg[0:len(Fht)]
        Fp_filt = p_reg[len(Fht):]
        Ff_filt = f_reg[len(Fht):]

        plt.plot(Fht, Fh_filt, '-', color="darkolivegreen", linewidth=Lwide*2, label='HIST')
        plt.plot(Fpt, Fp_filt, '-', color="brown", linewidth=Lwide*2, label='EVOMELT')
        plt.plot(Fpt, Ff_filt, '-', color="royalblue", linewidth=Lwide*2, label='FIXMELT')

    fsize = 8
    plt.xlim(np.array([np.min(Fht), np.max(Fpt)]))
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.tick_params(axis='both', labelsize=fsize)

    regname = Fh.transectNames[transect].values.item()
    plt.title(regname,fontsize=fsize)
    plt.xlabel('Time (a)', fontsize=fsize)
    plt.ylabel('Transport (Sv)', fontsize=fsize)
    plt.grid(True)
    plt.legend()
    plt.legend(fontsize=6)

    if opt_save == 1:
        plt.savefig(f'/Users/irenavankova/Work/data_sim/FISMF/Transport/Transp_{regname.replace(" ", "_")}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()
