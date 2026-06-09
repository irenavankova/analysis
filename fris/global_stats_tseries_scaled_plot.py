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

vars_global = ['CFLNumberGlobal', 'kineticEnergyCellMax', 'config_dt']
#ylabels_global = ['U from Global CFL (m/s)', 'EKE CellMax (m^2/s^-2)', 'Time step (s)']
#ylabels_global = [
#    r'U from Global CFL ($\mathrm{m\cdot s^{-1}}$)',
#    r'EKE CellMax ($\mathrm{m^2\cdot s^{-2}}$)',
#    r'Time step ($\mathrm{s}$)'
#]
ylabels_global = [
    r'Global CFL',
    r'EKE CellMax ($\mathrm{m^2\cdot s^{-2}}$)',
    r'Time step ($\mathrm{s}$)'
]
# Reference constant for CFL scaling: 60 sec / 1 km
dtdx_ref = 1000.0

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

print(f"\nGenerating 3-panel vertical column global stats plot...")

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(6, 8), dpi=150, sharex=True)

legend_handles = []
legend_labels = []

# Loop over files and extract pre-computed days
for cfg in valid_files:
    #with xr.open_dataset(cfg['path']) as ds:
    with xr.open_dataset(cfg['path'], decode_times=False) as ds:

        # Pull time axis directly from dataset coordinate/variable
        time_coords = ds['daysSinceStartOfSim'].values

        for idx, var_name in enumerate(vars_global):
            if var_name not in ds.data_vars:
                continue

            ax = axes[idx]
            plot_values = ds[var_name].values

            # Apply scale modifications to the CFL panel
            if var_name == 'CFLNumberGlobal':
                meshres = cfg['meshres']
                dt_array = ds['config_dt'].values
                #print(dt_array)
                #plot_values = (plot_values * meshres / dt_array) * dtdx_ref

            line, = ax.plot(time_coords, plot_values,
                            label=cfg['label'],
                            linewidth=1.0,
                            linestyle='-',
                            alpha=0.85,
                            color=cfg['color'])

            if idx == 0 and cfg['label'] not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(cfg['label'])

# Formatting panels
for idx, var_name in enumerate(vars_global):
    ax = axes[idx]
    ax.set_ylabel(ylabels_global[idx], fontsize=10)

    if idx == len(vars_global) - 1:
        ax.set_xlabel("Days Since Beginning of Simulation", fontsize=10)

    try:
        ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    except AttributeError:
        pass

if legend_handles:
    fig.legend(legend_handles, legend_labels, loc='upper center',
               bbox_to_anchor=(0.5, 0.95), ncol=len(legend_handles), frameon=True, fontsize=10)

plt.tight_layout(rect=[0, 0, 1, 0.91])

save_path = os.path.join(output_dir, "global_stats.png")
if opt_save == 1:
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Plot saved successfully to: {save_path}")
else:
    fig.set_dpi(95)
    plt.show()

plt.close(fig)