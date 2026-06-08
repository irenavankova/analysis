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

# Centralized line and time property definitions
Lwides1 = 1.0
Lwides6 = 1.0

lstyl1 = '--'
lstyl6 = '-'

S6_start_yr = 0

# Target dates for vertical lines (defined as tuple: (year, month))
vline_targets = [
    (2, 11),
    (0, 12),
    (4, 6)
]

vline_colors = ['yellowgreen', 'darkgray', 'darkgray']

# Filenames updated to match the global stats naming pattern
files_config = [
    {
        'filename': 'global_stats_tseries_F8_Spin1p1.nc',
        'label': 'F8 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'orange',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F8_Spin6p1.nc',
        'label': 'F8 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'brown',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'global_stats_tseries_F4_Spin1p1.nc',
        'label': 'F4 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'lightskyblue',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F4_Spin6p1.nc',
        'label': 'F4 (Spin6)',
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
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'forestgreen',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'global_stats_tseries_F2_Spin1p1.nc',
        'label': 'F2 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F2_Spin1p2.nc',
        'label': 'F2 (Spin1p2)',
        'start_year': 2,
        'start_month': 11,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F1_Spin1p1.nc',
        'label': 'F1 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F1_Spin1p2.nc',
        'label': 'F1 (Spin1p2)',
        'start_year': 0,
        'start_month': 12,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F1_Spin1p3.nc',
        'label': 'F1 (Spin1p3)',
        'start_year': 4,
        'start_month': 6,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'global_stats_tseries_F1_Spin6p1.nc',
        'label': 'F1 (Spin6)',
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

print(f"\nGenerating 3-panel global stats plot...")

# Expanded figsize from (15, 5) to (20, 5) to neatly fit the 3rd panel
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 5), dpi=150, sharex=True)

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

        # Plot each of the three global variables in their respective subplot
        for idx, var_name in enumerate(vars_global):
            if var_name not in ds.data_vars:
                continue

            ax = axes[idx]
            plot_values = ds[var_name].values

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

# Formatting and beautifying the subplots
for idx, var_name in enumerate(vars_global):
    ax = axes[idx]
    ax.set_title(var_name, fontsize=12, fontweight='bold', loc='left')
    ax.set_xlabel("Adjusted Model Year", fontsize=10)
    ax.set_ylabel(var_name, fontsize=10)

    # Only apply scalar formatting if the axis is actually numeric
    try:
        ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    except AttributeError:
        print(f"--> Note: Skipping scalar tick formatting for '{var_name}' (it contains categorical/string data).")

    # --- Draw the vertical marker lines ---
    for x_pos, line_color in zip(vline_positions, vline_colors):
        ax.axvline(x=x_pos, color=line_color, linestyle=':', linewidth=1.2, alpha=0.7, zorder=1)

fig.suptitle("Global Simulation Statistics Comparison", fontsize=16, fontweight='normal', y=1.02)

# Place the centralized legend across the top
if legend_handles:
    fig.legend(legend_handles, legend_labels, loc='upper center',
               bbox_to_anchor=(0.5, 0.96), ncol=len(legend_handles), frameon=True, fontsize=11)

plt.tight_layout(rect=[0, 0, 1, 0.90])

save_path = os.path.join(output_dir, "global_stats.png")
if opt_save == 1:
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Plot saved successfully to: {save_path}")
else:
    fig.set_dpi(95)
    plt.show()

plt.close(fig)