#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt


ff = 1
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

fname_lifw_hist = 'lifw_reg_hist_ave_tseries'
fname_lifw_pismf = 'lifw_reg_ave_tseries'
fname_hist = 'utsq_reg_hist_ave_tseries'
fname_pismf = 'utsq_reg_pismf_ave_tseries'
fname_fismf = 'utsq_reg_fismf_ave_tseries'

#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_lifw_hist = xr.open_dataset(f'{fpath}{fname_lifw_hist}.nc')
ds_lifw_pismf = xr.open_dataset(f'{fpath}{fname_lifw_pismf}.nc')
ds_hist = xr.open_dataset(f'{fpath}{fname_hist}.nc')
ds_pismf = xr.open_dataset(f'{fpath}{fname_pismf}.nc')
ds_fismf = xr.open_dataset(f'{fpath}{fname_fismf}.nc')

# Extract the DataArray (assuming variable is 'lifw')
#lifw = ds['lifw']  # dims: (Time, region)
lifw_h = ds_lifw_hist['lifw']  # dims: (Time, region)
n_time = len(ds_lifw_hist['Time'])
time_h = np.arange(0,n_time)/12 + 1950

lifw_p = ds_lifw_pismf['lifw']  # dims: (Time, region)
n_time = len(ds_lifw_pismf['Time'])
time_ssp = np.arange(0,n_time)/12 + 2015

utsq_h = ds_hist['utsq']  # dims: (Time, region)
utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f = ds_fismf['utsq']  # dims: (Time, region)

fs = 1/(time_ssp[1]-time_ssp[0])
fc = 1/((time_ssp[1]-time_ssp[0])*36)

# Plot time series for each region
for region in lifw_p.region.values:
    if region == "Antarctica":
        lifw_reg = lifw_p.sel(region=region)
        utp_reg = utsq_p.sel(region=region)
        kAnt = np.mean(lifw_reg.values / utp_reg.values)
        print(kAnt)

for region in lifw_p.region.values:
    plt.figure(figsize=(12, 6))

    lh_reg = lifw_h.sel(region=region)
    lp_reg = lifw_p.sel(region=region)
    uth_reg = utsq_h.sel(region=region)
    utp_reg = utsq_p.sel(region=region)
    utf_reg = utsq_f.sel(region=region)

    l_reg = np.append(lh_reg, lp_reg, axis=0)
    p_reg = np.append(uth_reg, utp_reg, axis=0)
    f_reg = np.append(uth_reg, utf_reg, axis=0)

    k = np.mean(lp_reg.values / utp_reg.values)

    print(k)
    Lwide = 1.25

    plt.plot(time_h, lh_reg, '-', color="black", linewidth=Lwide)
    plt.plot(time_h, uth_reg * kAnt, ':', color="yellowgreen", linewidth=Lwide)

    plt.plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide)
    plt.plot(time_ssp, utp_reg * kAnt, ':', color="darkorange", linewidth=Lwide)

    plt.plot(time_ssp, utf_reg * kAnt, ':', color="lightskyblue", linewidth=Lwide)

    if ff == 1:
        l_reg = iv_filt.butter_filter(l_reg, fc, fs, 'low')
        p_reg = iv_filt.butter_filter(p_reg, fc, fs, 'low')
        f_reg = iv_filt.butter_filter(f_reg, fc, fs, 'low')

        lh_reg = l_reg[0:len(time_h)]
        lp_reg = l_reg[len(time_h):]

        uth_reg = p_reg[0:len(time_h)]
        utp_reg = p_reg[len(time_h):]

        utf_reg = f_reg[len(time_h):]

        plt.plot(time_h, lh_reg, '-', color="black", linewidth=Lwide, label='lifw_p')
        plt.plot(time_h, uth_reg * kAnt, '--', color="darkolivegreen", linewidth=Lwide, label='ut_p')

        plt.plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide, label='lifw_p')
        plt.plot(time_ssp, utp_reg * kAnt, '--', color="brown", linewidth=Lwide, label='ut_p')

        plt.plot(time_ssp, utf_reg * kAnt, '-', color="royalblue", linewidth=Lwide, label='ut_f')

    plt.title(region)
    plt.xlabel('Time')
    plt.grid(True)
    plt.legend()
    plt.show()
