#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import mode
import datetime
import gsw

# %% OBSERVATIONS
P = np.array([1125, 1000, 910, 820])
z = gsw.z_from_p(P, -78)

# Load data file
A = np.loadtxt('/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/DOVuFRIS/Fris_Apres/code/ocean_data/Site3_ocean_data/Moorings/ifixread.all')

# %% Processing and Filtering
time = A[:, 0]
ii = np.where((time > 26) & (time < 1015))[0]

dt = mode(np.diff(time), keepdims=True).mode[0]
tt = np.arange(time[ii[0]], time[ii[-1]] + dt, dt)

# Convert relative time 'tt' to datetime objects (1996 base year)
base_date = datetime.datetime(1996, 1, 1)
td = np.array([base_date + datetime.timedelta(days=float(t)) for t in tt])

fc = 1 / 60  # Low-pass filter cutoff
fs = 1 / dt  # Sampling frequency
b, a = butter(4, fc / (fs / 2), btype='low')

T_obs = np.zeros((4, len(tt)))

# --- Plot 1: Time Series Lines ---
plt.figure(figsize=(10, 4))
for j in range(4):
    col_idx = (j + 1) * 3 - 1
    y_med = medfilt(A[ii, col_idx], kernel_size=5)
    _, unique_indices = np.unique(time[ii], return_index=True)

    f_interp = interp1d(time[ii[unique_indices]], y_med[unique_indices], kind='linear', fill_value="extrapolate")
    y_interp = f_interp(tt)

    yy = filtfilt(b, a, y_interp)
    T_obs[j, :] = yy

    plt.plot(td, yy)

plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%m'))
plt.gcf().autofmt_xdate()
plt.grid(True)
plt.title('Observed Temperature Time Series')

# --- Plot 2: Observation Contour with hlines ---
hlines = np.arange(-2.45, -1.91 + 0.01, 0.03)

fig, ax = plt.subplots(figsize=(10, 6))

# 'bone_r' is the equivalent to MATLAB's flipud(bone)
cf = ax.contourf(td, z, T_obs, levels=hlines, cmap='hot_r')
cl = ax.contour(td, z, T_obs, levels=hlines, colors='k', linewidths=0.5)

# FIX 1: set_clim belongs to the contour object (cf), not the axes (ax)
cf.set_clim(min(hlines), max(hlines))

ax.grid(True, color='white')
ax.set_facecolor('#EAEAEA')

# Set specific X ticks (Jan 16th for years 1996 to 1999)
tick_years = [datetime.datetime(y, 1, 16, 12, 0, 0) for y in range(1996, 2000)]
ax.set_xticks(tick_years)
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y'))

ax.set_ylabel('z (m)')

# FIX 2: Added 'r' for a raw string to fix the LaTeX "\c" warning
ax.set_title(r'Obs: Site 3, in situ temperature ($^\circ$C)', fontweight='normal', fontsize=8)
ax.tick_params(labelsize=8)

# Colorbar with explicit hlines ticks
cb = fig.colorbar(cf, ax=ax, ticks=hlines)
cb.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.show()