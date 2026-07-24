#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import gsw

# %% OBSERVATIONS TARGET DEPTHS
P = np.array([1125, 1000, 910, 820])
zobs = gsw.z_from_p(P, -78)

# ==========================================
# 1. DEFINE MODEL CONFIGURATIONS
# ==========================================
MODEL_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/obs"
runs = ['F1', 'F2', 'F4', 'F8']

MODEL_THICK_VAR = "timeMonthly_avg_layerThickness"
MODEL_SSH_VAR = "timeMonthly_avg_ssh"
# Switched from temperature to salinity variable
SALT_VAR = "timeMonthly_avg_activeTracers_salinity"

# ==========================================
# 2. SET UP 2x2 SUBPLOT MATRIX
# ==========================================
fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)
axes = axes.flatten()  # 1D array of 4 elements

# Use a consistent color palette for the 4 target depths across all subplots
colors = plt.cm.tab10(np.linspace(0, 1, len(zobs)))

# --- Plot Model Configurations ---
for run_idx, run in enumerate(runs):
    ax = axes[run_idx]

    file_name = f'obs_tseries_{run}_Spin6p1.nc'
    file_path = os.path.join(MODEL_DIRECTORY, file_name)

    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f"File not found:\n{file_name}", ha='center', va='center', color='red')
        ax.set_title(f"Model Run: {run} (Missing)", fontweight='bold')
        continue

    ds = xr.open_dataset(file_path)

    # Extract Site 3 coordinates
    salt_site3 = ds[SALT_VAR].sel(site='Site3')
    thick_site3 = ds[MODEL_THICK_VAR].sel(site='Site3')
    ssh_site3 = ds[MODEL_SSH_VAR].sel(site='Site3')

    # Filter out empty vertical layers based on salinity data availability
    salt_filtered = salt_site3.dropna(dim='nVertLevels', how='all')
    valid_levels = salt_filtered['nVertLevels'].values

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

        # Generate custom label strings exclusively on the final plot pass for the global legend
        label_text = f"Target {zo:.0f}m (Lvl {lvl}: {actual_model_depth:.1f}m)" if run_idx == 3 else ""

        ax.plot(
            model_years,
            salt_filtered.sel(nVertLevels=lvl).values,
            color=colors[i],
            linewidth=1.8,
            label=label_text
        )

    ax.set_title(f'Model Run: {run}', fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ds.close()

# ==========================================
# 3. GRID-WIDE LABELS & STYLING
# ==========================================
# Set X labels for the bottom row subplots (indices 2, 3)
for idx in [2, 3]:
    axes[idx].set_xlabel('Model Year (Starting at 0)', fontsize=11)

# Set Y labels for the leftmost column subplots (indices 0, 2)
for idx in [0, 2]:
    axes[idx].set_ylabel('Salinity (g/kg or psu)', fontsize=11)

# Grab the line handle information from the completed loop to form a unified bottom legend
handles, labels = axes[3].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_transform=fig.transFigure, bbox_to_anchor=(0.5, 0.01), ncol=2, frameon=True)

fig.suptitle('Site 3: Model Salinity Inter-comparison (F1, F2, F4, F8)', fontsize=16, fontweight='bold', y=0.96)

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
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, "comp_obs_MoorSite3_S.png"),
                dpi=OUTPUT_DPI,
                bbox_inches='tight')

plt.show()