#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================
data_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/global_stats/par4dt'
opt_save = 1

vars_global = ['CFLNumberGlobal', 'kineticEnergyCellMax']

ylabels_global = [
    r'speed of max CFL ($\mathrm{m\cdot s^{-1}}$)',
    r'max cell kinetic energy ($\mathrm{m^2\cdot s^{-2}}$)',
    r'Time step ($\mathrm{s}$)'
]
# Reference constant for CFL scaling: 60 sec / 1 km
dtdx_ref = 1000.0
km2m = 1000
# Filenames config
files_config = [
    {'filename': 'global_stats_tseries_F8_Spin6p1.nc', 'label': 'F8', 'meshres': 8.0, 'color': 'brown'},
    {'filename': 'global_stats_tseries_F4_Spin6p1.nc', 'label': 'F4', 'meshres': 4.0, 'color': 'dodgerblue'},
    {'filename': 'global_stats_tseries_F2_Spin6p1.nc', 'label': 'F2', 'meshres': 2.0, 'color': 'orange'},
    {'filename': 'global_stats_tseries_F1_Spin6p1.nc', 'label': 'F1', 'meshres': 1.0, 'color': 'black'},
]

output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/global_stats/'
os.makedirs(output_dir, exist_ok=True)

# ==============================================================================

print("Validating explicitly configured files...")
valid_files = []
for cfg in files_config:
    full_path = os.path.join(data_dir, cfg['filename'])
    if os.path.exists(full_path):
        cfg['path'] = full_path
        valid_files.append(cfg)
        print(f"--> Confirmed: {cfg['filename']} | Label: {cfg['label']}")
    else:
        print(f"Warning: File not found, skipping: {full_path}")

if not valid_files:
    raise FileNotFoundError("Error: None of the explicitly specified NetCDF files were found.")

# DYNAMIC ROW AND HEIGHT CALCULATION
num_vars = len(vars_global)
height_per_panel = 2.66
fig_height = max(3.0, num_vars * height_per_panel)

print(f"\nGenerating {num_vars}-panel vertical column global stats plot (Height: {fig_height:.2f} inches)...")

fig, axes = plt.subplots(nrows=num_vars, ncols=1, figsize=(6, fig_height), dpi=150, sharex=True, squeeze=False)
axes = axes.flatten()

legend_handles = []
legend_labels = []

# Tracking global min/max years to set tight x-axis limits later
all_years_min = []
all_years_max = []

# Loop over files and extract pre-computed days
for cfg in valid_files:
    with xr.open_dataset(cfg['path'], decode_times=False) as ds:

        # Convert days axis to years axis directly
        time_coords_days = ds['daysSinceStartOfSim'].values
        time_coords_years = time_coords_days / 365.25

        all_years_min.append(np.nanmin(time_coords_years))
        all_years_max.append(np.nanmax(time_coords_years))

        for idx, var_name in enumerate(vars_global):
            if var_name not in ds.data_vars:
                continue

            ax = axes[idx]
            plot_values = ds[var_name].values

            # Apply scale modifications to the CFL panel
            if var_name == 'CFLNumberGlobal':
                meshres = cfg['meshres'] * km2m
                dt_array = ds['config_dt'].values
                #print(meshres)
                #print(dt_array)
                plot_values = (plot_values * meshres / dt_array)

            line, = ax.plot(time_coords_years, plot_values,
                            label=cfg['label'],
                            linewidth=1.0,
                            linestyle='-',
                            alpha=0.85,
                            color=cfg['color'])

            if idx == 0 and cfg['label'] not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(cfg['label'])

# Global tight bounds calculation
x_min = min(all_years_min) if all_years_min else 0
x_max = max(all_years_max) if all_years_max else 1

# Formatting panels
for idx, var_name in enumerate(vars_global):
    ax = axes[idx]
    ax.set_ylabel(ylabels_global[idx], fontsize=10)

    # Force tight limits on the x-axis for every subplot
    ax.set_xlim(x_min, x_max)

    # Add subplot labels (a), (b), (c)... to the upper left corner
    panel_label = f"{chr(97 + idx)})"  # Converts 0 to (a), 1 to (b), etc.
    ax.text(0.01, 0.95, panel_label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')

    if idx == len(vars_global) - 1:
        ax.set_xlabel("RUN years", fontsize=10)

    try:
        ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    except AttributeError:
        pass

# Adjust subplots tight boundary constraints first
plt.tight_layout()

# Draw the legend strictly in the upper right corner of the figure canvas
if legend_handles:
    fig.legend(
        legend_handles,
        legend_labels,
        loc='upper right',
        bbox_to_anchor=(0.98, 0.98),
        ncol=1,  # 1 column layout stack (rows)
        frameon=True,
        fontsize=8  # 2 pts smaller than the axis labels (10 pts)
    )

save_path = os.path.join(output_dir, "global_cfl_u.png")
if opt_save == 1:
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Plot saved successfully to: {save_path}")
else:
    fig.set_dpi(95)
    plt.show()

plt.close(fig)