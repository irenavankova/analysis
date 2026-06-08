#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt

# ==============================================================================
# USER CONFIGURATION: Explicitly specify your files, labels, and start targets
# ==============================================================================
data_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/global_stats/par4dt'
opt_save = 0

# Variables to plot in the 3 subplots
vars_global = ['CFLNumberGlobal', 'kineticEnergyCellMax', 'config_dt']
ylabels_global = ['Scaled Global CFL', 'EKE CellMax (m^2 s^-2)', 'dt (s)']

# Reference constant for CFL scaling: 60 sec / 1 km
dtdx_ref = 60.0 / 1.0

# Centralized line and time property definitions
Lwides6 = 1.0
lstyl6 = '-'
S6_start_yr = 0

# Target dates for vertical lines (defined as tuple: (year, month))
vline_targets = [
    (2, 11),
    (0, 12),
    (4, 6)
]

vline_colors = ['yellowgreen', 'darkgray', 'darkgray']

# Filenames config filtered to ONLY contain Spin6 runs
files_config = [
    {
        'filename': 'global_stats_tseries_F8_Spin6p1.nc',
        'label': 'F8 (Spin6)',
        'meshres': 8.0,
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'brown',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'global_stats_tseries_F4_Spin6p1.nc',
        'label': 'F4 (Spin6)',
        'meshres': 4.0,
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'royalblue',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'global_stats_tseries_F2_Spin6p1.nc',
        'label': 'F2 (Spin6)',
        'meshres': 2.0,
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'forestgreen',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'global_stats_tseries_F1_Spin6p1.nc',
        'label': 'F1 (Spin6)',
        'meshres': 1.0,
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'black',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
]

output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/global_stats/'
os.makedirs(output_dir, exist_ok=True)

# ==============================================================================

# Pre-calculate the decimal positions for your vertical lines
vline_positions = []
for yr, mo in vline_targets:
    dec_yr = yr + ((mo - 1) / 12.0)
    vline_positions.append(dec_yr)

print("Validating explicitly configured files...")
valid_files = []
for cfg in files_config:
    full_path = os.path.join(data_dir, cfg['filename'])
    if os.path.exists(full_path):
        cfg['path'] = full_path
        valid_files.append(cfg)
        print(f"--> Confirmed: {cfg['filename']} | Label: {cfg['label']}")
    else:
        print(f"Warning: Explicitly specified file not found, skipping: {full_path}")

if not valid_files:
    raise FileNotFoundError("Error: None of the explicitly specified NetCDF files were found.")

print(f"\nGenerating 3-panel vertical column global stats plot (Spin6 Only)...")

# 3 rows, 1 column layout
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), dpi=150, sharex=True)

legend_handles = []
legend_labels = []

# Loop over each valid netCDF file and plot data on all axes
for cfg in valid_files:
    with xr.open_dataset(cfg['path']) as ds:

        # Get the length of your time dimension from the dataset
        num_timesteps = ds.sizes['Time']
        start_year = cfg['start_year']
        start_month = cfg['start_month']

        time_coords = []
        for i in range(num_timesteps):
            total_months = (start_month - 1) + i
            current_year = start_year + (total_months // 12)
            current_month_0indexed = total_months % 12
            decimal_year = current_year + (current_month_0indexed / 12.0)
            time_coords.append(decimal_year)

        # Plot each of the three global variables in their respective subplot row
        for idx, var_name in enumerate(vars_global):
            if var_name not in ds.data_vars:
                continue

            ax = axes[idx]
            plot_values = ds[var_name].values

            # Apply scale modifications to the CFL panel
            if var_name == 'CFLNumberGlobal':
                meshres = cfg['meshres']
                dt_array = ds['config_dt'].values
                plot_values = (plot_values * meshres / dt_array) * dtdx_ref

            line, = ax.plot(time_coords, plot_values,
                            label=cfg['label'],
                            linewidth=cfg.get('linewidth', 1.5),
                            linestyle=cfg.get('linestyle', '-'),
                            alpha=0.85,
                            color=cfg.get('color', None))

            # Fill legend tracking from the first subplot to avoid duplicates
            if idx == 0 and cfg.get('in_legend', True) and cfg['label'] not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(cfg['label'])

# Formatting and beautifying the stacked subplots
for idx, var_name in enumerate(vars_global):
    ax = axes[idx]

    #display_title = f"{var_name} (Scaled)" if var_name == 'CFLNumberGlobal' else var_name

    #ax.set_title(display_title, fontsize=12, fontweight='bold', loc='left')
    ax.set_ylabel(ylabels_global, fontsize=10)

    # Only label the bottom panel's x-axis to prevent overlapping layout clutter
    if idx == len(vars_global) - 1:
        ax.set_xlabel("Adjusted Model Year", fontsize=10)

    # Only apply scalar formatting if the axis is actually numeric
    try:
        ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    except AttributeError:
        pass

    # --- Draw the vertical marker lines ---
    #for x_pos, line_color in zip(vline_positions, vline_colors):
    #    ax.axvline(x=x_pos, color=line_color, linestyle=':', linewidth=1.2, alpha=0.7, zorder=1)

fig.suptitle("Global Simulation Statistics Comparison (Spin6 runs)", fontsize=16, fontweight='normal', y=0.98)

# Centralized legend adjusted across the top frame (since there are fewer entries, ncol=len entries fits well)
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