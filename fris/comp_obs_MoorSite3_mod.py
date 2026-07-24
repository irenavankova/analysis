#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import datetime
import gsw
from scipy.interpolate import interp1d
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import mode

# ==========================================
# 1. PROCESS OBSERVATIONS
# ==========================================
P = np.array([1125, 1000, 910, 820])
zobs = gsw.z_from_p(P, -78)

# Load data file
obs_file = '/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/DOVuFRIS/Fris_Apres/code/ocean_data/Site3_ocean_data/Moorings/ifixread.all'
A = np.loadtxt(obs_file)

time = A[:, 0]
ii = np.where((time > 26) & (time < 1015))[0]

dt = mode(np.diff(time), keepdims=True).mode[0]
tt = np.arange(time[ii[0]], time[ii[-1]] + dt, dt)

# Filter parameters
fc = 1 / 60  # Low-pass filter cutoff
fs = 1 / dt  # Sampling frequency
b, a = butter(4, fc / (fs / 2), btype='low')

T_obs = np.zeros((4, len(tt)))

for j in range(4):
    col_idx = (j + 1) * 3 - 1
    y_med = medfilt(A[ii, col_idx], kernel_size=5)
    _, unique_indices = np.unique(time[ii], return_index=True)

    f_interp = interp1d(time[ii[unique_indices]], y_med[unique_indices], kind='linear', fill_value="extrapolate")
    y_interp = f_interp(tt)
    T_obs[j, :] = filtfilt(b, a, y_interp)

# Shift observation time to start at Year 0 (Decimal Years)
obs_days_since_start = tt - tt[0]
obs_years = obs_days_since_start / 365.0

# ==========================================
# 2. DEFINE MODEL CONFIGURATIONS
# ==========================================
MODEL_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/obs"
runs = ['F1', 'F2', 'F4', 'F8']

MODEL_THICK_VAR = "timeMonthly_avg_layerThickness"
MODEL_SSH_VAR = "timeMonthly_avg_ssh"
TEMP_VAR = "timeMonthly_avg_activeTracers_temperature"

# ==========================================
# 3. SET UP 2x3 SUBPLOT MATRIX
# ==========================================
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()  # 1D array of 6 elements

# Use a consistent color palette for the 4 target depths across all subplots
colors = plt.cm.tab10(np.linspace(0, 1, len(zobs)))

# --- Plot Observations in Subplot 0 ---
ax_obs = axes[0]
for i in range(4):
    label_text = f"Target {zobs[i]:.0f}m"
    ax_obs.plot(obs_years, T_obs[i, :], color=colors[i], linewidth=1.8, label=label_text)

ax_obs.set_title('Observations (MoorSite3)', fontsize=12, fontweight='bold')
ax_obs.grid(True, linestyle='--', alpha=0.5)
ax_obs.set_ylabel('Potential Temperature (°C)', fontsize=11)

# --- Plot Model Configurations in Subplots 1 through 4 ---
for run_idx, run in enumerate(runs):
    # Shift index by 1 because observations occupy the first subplot index [0]
    ax = axes[run_idx + 1]

    file_name = f'obs_tseries_{run}_Spin6p1.nc'
    file_path = os.path.join(MODEL_DIRECTORY, file_name)

    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f"File not found:\n{file_name}", ha='center', va='center', color='red')
        ax.set_title(f"Model Run: {run} (Missing)", fontweight='bold')
        continue

    ds = xr.open_dataset(file_path)

    # Extract Site 3 coordinates
    temp_site3 = ds[TEMP_VAR].sel(site='Site3')
    thick_site3 = ds[MODEL_THICK_VAR].sel(site='Site3')
    ssh_site3 = ds[MODEL_SSH_VAR].sel(site='Site3')

    # Filter out empty vertical layers
    temp_filtered = temp_site3.dropna(dim='nVertLevels', how='all')
    valid_levels = temp_filtered['nVertLevels'].values

    # Calculate static time-averaged z profile structure
    mean_ssh = float(ssh_site3.mean(dim='Time').values)
    mean_thickness = thick_site3.mean(dim='Time').values

    mean_z = np.zeros(len(mean_thickness))
    current_top = mean_ssh
    for k in range(len(mean_thickness)):
        mean_z[k] = current_top - 0.5 * mean_thickness[k]
        current_top -= mean_thickness[k]

    mean_z_filtered = mean_z[valid_levels]

    # Convert model time to decimal Model Years (Year 0 start)
    time_days = ds['Time'].values
    timedeltas = time_days - time_days[0]
    days_since_start = np.array([td.days + (td.seconds / 86400.0) for td in timedeltas])
    model_years = days_since_start / 365.0

    # Match and plot closest model layer for each target zobs
    for i, zo in enumerate(zobs):
        idx = np.abs(mean_z_filtered - zo).argmin()
        lvl = valid_levels[idx]
        actual_model_depth = mean_z_filtered[idx]

        # Pull distinct descriptive metadata only for the shared global legend loop
        label_text = f"Target {zo:.0f}m (Lvl {lvl}: {actual_model_depth:.1f}m)" if run_idx == 3 else ""

        ax.plot(
            model_years,
            temp_filtered.sel(nVertLevels=lvl).values,
            color=colors[i],
            linewidth=1.8,
            label=label_text
        )

    ax.set_title(f'Model Run: {run}', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ds.close()

# --- Hide the extra 6th subplot pane ---
axes[5].axis('off')

# ==========================================
# 4. GRID-WIDE LABELS & STYLING
# ==========================================
# Set X labels for the bottom row subplots (indices 3, 4)
for idx in [3, 4]:
    axes[idx].set_xlabel('Elapsed Year (Starting at 0)', fontsize=11)

# Set Y labels for the leftmost column subplots (indices 0, 3)
for idx in [0, 3]:
    axes[idx].set_ylabel('Potential Temperature (°C)', fontsize=11)

# Draw a single clean layout legend at the bottom of the entire matrix
# Using labels from the observation axes ensures basic target tracking is clear
handles, labels = ax_obs.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_transform=fig.transFigure, bbox_to_anchor=(0.5, 0.01), ncol=4,
           frameon=True)

fig.suptitle('Site 3: Observation vs. Model Temperature Inter-comparison (F1, F2, F4, F8)', fontsize=16,
             fontweight='bold', y=0.96)

plt.tight_layout(rect=[0, 0.08, 1, 0.95])

# ==========================================
# 5. SAVE AND SHOW
# ==========================================
SAVE_FIGURE = True
OUTPUT_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_plots/compare_obs"
OUTPUT_DPI = 300

if SAVE_FIGURE:
    # Automatically create the directory if it doesn't exist
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

    # Save the figure safely
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, "comp_obs_MoorSite3_T.png"),
                dpi=OUTPUT_DPI,
                bbox_inches='tight')

plt.show()