#!/usr/bin/env python3
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import iv_filt
from scipy.optimize import curve_fit


ff = 1
fpath = '/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/post_derived/'

fname_lifw_hist = 'lifw_reg_hist_ave_tseries'
fname_lifw_pismf = 'lifw_reg_ave_tseries'
fname_lifw_fismf = 'lifw_reg_fismf_ave_tseries'

fname_hist = 'utsq_reg_hist_ave_tseries'
fname_pismf = 'utsq_reg_pismf_ave_tseries'
fname_fismf = 'utsq_reg_fismf_ave_tseries'

#ds = xr.open_dataset(f'{fpath}{fname}.nc')
ds_lifw_hist = xr.open_dataset(f'{fpath}{fname_lifw_hist}.nc')
ds_lifw_pismf = xr.open_dataset(f'{fpath}{fname_lifw_pismf}.nc')
ds_lifw_fismf = xr.open_dataset(f'{fpath}{fname_lifw_fismf}.nc')

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

lifw_f = ds_lifw_fismf['lifw']  # dims: (Time, region)


utsq_h = ds_hist['utsq']  # dims: (Time, region)
utsq_p = ds_pismf['utsq']  # dims: (Time, region)
utsq_f = ds_fismf['utsq']  # dims: (Time, region)

fs = 1/(time_ssp[1]-time_ssp[0])
fc = 1/((time_ssp[1]-time_ssp[0])*36)

def f_exp(x, a,b,c):
    return a * np.exp(x/b) + c

def f_n_lin(x, a, b):
    return a * x + b


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
    lf_reg = lifw_f.sel(region=region)

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

        plt.plot(time_h, lh_reg, '-', color="black", linewidth=Lwide+2)
        plt.plot(time_h, uth_reg * kAnt, '--', color="darkolivegreen", linewidth=Lwide+2, label='HIST')

        plt.plot(time_ssp, lp_reg, '-', color="black", linewidth=Lwide+2)
        plt.plot(time_ssp, utp_reg * kAnt, '--', color="brown", linewidth=Lwide+2, label='EVOMELT')

        plt.plot(time_ssp, utf_reg * kAnt, '-', color="royalblue", linewidth=Lwide+2, label='FIXMELT')

        popt, pcov = curve_fit(f_n_lin, time_ssp-time_ssp[0], utf_reg * kAnt)
        plt.plot(time_ssp, f_n_lin(time_ssp-time_ssp[0], *popt), '-', color="magenta", linewidth=2)

        popt, pcov = curve_fit(f_exp, time_ssp - time_ssp[0], utf_reg * kAnt, p0=[1, 100, 1200])
        plt.plot(time_ssp, f_exp(time_ssp - time_ssp[0], *popt), '--', color="cyan", linewidth=2)

        offset = time_ssp[0]
        offset = 2040
        i1 = 12*25
        popt, pcov = curve_fit(f_n_lin, time_ssp[i1:] - offset, utp_reg[i1:] * kAnt)
        plt.plot(time_ssp[i1:], f_n_lin(time_ssp[i1:] - offset, *popt), '-', color="magenta", linewidth=2)

        popt, pcov = curve_fit(f_exp, time_ssp[i1:] - offset, utp_reg[i1:] * kAnt, p0=[1, 100, 1200])
        plt.plot(time_ssp[i1:], f_exp(time_ssp[i1:] - offset, *popt), '--', color="cyan", linewidth=2)

        fflux = np.sum(utf_reg * kAnt)/12 - np.sum(lf_reg.values)/12
        pflux = np.sum(utp_reg * kAnt)/12 - np.sum(lf_reg.values)/12

        print(fflux)

        # Calculate the center of the plot
        x_center = (plt.xlim()[0] + plt.xlim()[1]) / 2
        y_center = (plt.ylim()[0] + plt.ylim()[1]) / 2

        # Add text to the center of the plot
        #plt.text(x_center, y_center, f'Ffulx = {fflux} GT', ha='center', va='center', fontsize=12, color='red')
        #plt.text(x_center, y_center+y_center/8, f'Pflux = {pflux} GT', ha='center', va='center', fontsize=12, color='red')
        #plt.text(x_center, y_center+y_center/4, f'diff = {pflux-fflux} GT', ha='center', va='center', fontsize=12, color='red')
        #plt.text(x_center, y_center+y_center/2, f'perc = {(pflux-fflux)/fflux*100} %', ha='center', va='center', fontsize=12, color='red')


    plt.plot(time_ssp, lf_reg, '-', color="gold", linewidth=Lwide+2)

    plt.title(region)
    plt.xlabel('Time')
    plt.grid(True)
    plt.legend()
    plt.show()
