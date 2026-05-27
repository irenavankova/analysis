#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt

# ==============================================================================
# USER CONFIGURATION: Explicitly specify your files, labels, and start targets
# ==============================================================================
# Directory where all your netCDF files are located ('.' means current directory)
data_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/bulk_seaice_tseries'
opt_save = 1

# Centralized line and time property definitions
Lwides1 = 2.0
Lwides6 = 1.5

lstyl1 = '--'
lstyl6 = '-'

S6_start_yr = 0

# Target dates for vertical lines (defined as tuple: (year, month))
vline_targets = [
    (2, 11),
    (0, 12),
    (4, 6)
]

# Note: Ensure these filenames match the sea-ice netCDF outputs you generated
# (e.g., bulk_seaice_tseries_...)
files_config = [
    {
        'filename': 'bulk_seaice_tseries_F8_Spin1.nc',
        'label': 'F8 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'orange',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F8_Spin6.nc',
        'label': 'F8 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'brown',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_seaice_tseries_F4_Spin1.nc',
        'label': 'F4 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'lightskyblue',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F4_Spin6.nc',
        'label': 'F4 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'royalblue',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_seaice_tseries_F2_Spin6.nc',
        'label': 'F2 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'forestgreen',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_seaice_tseries_F2_Spin1p1.nc',
        'label': 'F2 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F2_Spin1p2.nc',
        'label': 'F2 (Spin1p2)',
        'start_year': 2,
        'start_month': 11,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F1_Spin1p1.nc',
        'label': 'F1 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F1_Spin1p2.nc',
        'label': 'F1 (Spin1p2)',
        'start_year': 0,
        'start_month': 12,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F1_Spin1p3.nc',
        'label': 'F1 (Spin1p3)',
        'start_year': 4,
        'start_month': 6,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_seaice_tseries_F1_Spin6.nc',
        'label': 'F1 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'black',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
]

output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/seaice_tseries/compare_all/'
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

# Collect distinct variables and regions
plot_vars = set()
regions = []
for cfg in valid_files:
    with xr.open_dataset(cfg['path']) as ds:
        vars_in_file = [v for v in ds.data_vars if v not in ['Time', 'region']]
        plot_vars.update(vars_in_file)
        if not regions:
            regions = ds['region'].values.tolist()

plot_vars = sorted(list(plot_vars))
num_regions = len(regions)
print(f"\nProcessing {len(plot_vars)} variables across {num_regions} regions...")

ncols = 2
nrows = (num_regions + 1) // ncols

for var_name in plot_vars:
    print(f"Generating multi-panel comparison plot for: {var_name}...")

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 3.5 * nrows), dpi=150, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    legend_handles = []
    legend_labels = []

    for cfg in valid_files:
        with xr.open_dataset(cfg['path']) as ds:

            if var_name not in ds.data_vars:
                continue

            start_year = cfg['start_year']
            start_month = cfg['start_month']

            # 2. Get the length of your time dimension from the dataset
            num_timesteps = ds.sizes['Time']  # e.g., 11

            time_coords = []
            for i in range(num_timesteps):
                # Calculate total elapsed months from the target starting point
                # (Subtracting 1 shifts month from 1-12 down to a 0-11 index for clean math)
                total_months = (start_month - 1) + i

                # Use floor division and modulo to cleanly advance years and months
                current_year = start_year + (total_months // 12)
                current_month_0indexed = total_months % 12

                # Convert it back to your decimal year format
                decimal_year = current_year + (current_month_0indexed / 12.0)
                time_coords.append(decimal_year)

            data_slice = ds[var_name]

            for reg_idx, region in enumerate(regions):
                ax = axes_flat[reg_idx]

                if region in ds['region'].values:
                    reg_data = data_slice.sel(region=region)
                    plot_values = reg_data.values

                    # Unit conversions if necessary
                    # (Note: production/melting are scaled to m/yr directly in your analysis script)
                    if 'concentration' in var_name:
                        # Convert fractional area fraction (0-1) to percentage if preferred
                        # plot_values = plot_values * 100.0
                        pass

                    line, = ax.plot(time_coords, plot_values,
                                    label=cfg['label'],
                                    linewidth=cfg.get('linewidth', 1.5),
                                    linestyle=cfg.get('linestyle', '-'),
                                    alpha=0.85,
                                    color=cfg.get('color', None))

                    if cfg.get('in_legend', True) and cfg['label'] not in legend_labels:
                        legend_handles.append(line)
                        legend_labels.append(cfg['label'])

    # Step 5: Formatting, beautifying, and adding vertical lines
    for reg_idx, region in enumerate(regions):
        ax = axes_flat[reg_idx]
        ax.set_title(f"Region: {region}", fontsize=11, fontweight='normal', loc='left')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.ticklabel_format(style='plain', useOffset=False, axis='y')

        # --- Draw the vertical marker lines ---
        for x_pos in vline_positions:
            ax.axvline(x=x_pos, color='red', linestyle=':', linewidth=1.2, alpha=0.7, zorder=1)

        if reg_idx % ncols == 0:  # Only label leftmost y-axes
            if 'concentration' in var_name:
                ax.set_ylabel("Mean Ice Concentration (Fraction)", fontsize=10)
            elif 'production' in var_name:
                ax.set_ylabel("Mean Sea Ice Production ($m/year$)", fontsize=10)
            elif 'melting' in var_name:
                ax.set_ylabel("Mean Sea Ice Melting ($m/year$)", fontsize=10)
            elif 'integrated_ice_volume' in var_name:
                ax.set_ylabel("Integrated Ice Volume ($m^3$)", fontsize=10)
            elif 'thickness' in var_name:
                if 'max' in var_name:
                    ax.set_ylabel("Max Ice Thickness ($m$)", fontsize=10)
                else:
                    ax.set_ylabel("Mean Ice Thickness ($m$)", fontsize=10)
            else:
                ax.set_ylabel("Value", fontsize=10)

        if reg_idx >= (num_regions - ncols):
            ax.set_xlabel("Adjusted Model Year", fontsize=10)

    # Clean up leftover empty subplot blocks
    for residual_idx in range(num_regions, len(axes_flat)):
        fig.delaxes(axes_flat[residual_idx])

    fig.suptitle(f"Regional Sea Ice Comparison: {var_name}", fontsize=16, fontweight='normal', y=0.99)

    if legend_handles:
        fig.legend(legend_handles, legend_labels, loc='upper center',
                   bbox_to_anchor=(0.5, 0.96), ncol=len(legend_handles), frameon=True, fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    save_path = os.path.join(output_dir, f"compare_{var_name}.png")
    if opt_save == 1:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        fig.set_dpi(95)
        plt.show()
    plt.close(fig)

print(f"\nSuccess! All multi-panel sea-ice comparison plots completed.")