#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt


# ==============================================================================
# USER CONFIGURATION: Explicitly specify your files, labels, and start targets
# ==============================================================================
# Directory where all your netCDF files are located ('.' means current directory)
data_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/bulk_2D_tseries/'
opt_save = 1

# Centralized line and time property definitions
Lwides1 = 2.0
Lwides6 = 1.5

lstyl1 = '--'
lstyl6 = '-'

S6_start_yr = 0

files_config = [
    {
        'filename': 'bulk_tseries_2D_1by1_F8_Spin1.nc',
        'label': 'F8 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'orange',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F8_Spin6.nc',
        'label': 'F8 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'brown',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F4_Spin1.nc',
        'label': 'F4 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'lightskyblue',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F4_Spin6.nc',
        'label': 'F4 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'royalblue',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F2_Spin6.nc',
        'label': 'F2 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'forestgreen',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F2_Spin1p1.nc',
        'label': 'F2 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
{
        'filename': 'bulk_tseries_2D_1by1_F2_Spin1p2.nc',
        'label': 'F2 (Spin1p2)',
        'start_year': 2,
        'start_month': 11,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p1.nc',
        'label': 'F1 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p2.nc',
        'label': 'F1 (Spin1p2)',
        'start_year': 0,
        'start_month': 12,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p3.nc',
        'label': 'F1 (Spin1p3)',
        'start_year': 4,
        'start_month': 6,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin6.nc',
        'label': 'F1 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'black',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True
    },
]

output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/melt_tseries/compare_all/'
os.makedirs(output_dir, exist_ok=True)

rho_fw = 1000.
sec_per_day = 86400.
sec_per_year = sec_per_day * 365.

# ==============================================================================

print("Validating explicitly configured files...")
valid_files = []
for cfg in files_config:
    full_path = os.path.join(data_dir, cfg['filename'])
    if os.path.exists(full_path):
        cfg['path'] = full_path
        valid_files.append(cfg)
        print(
            f"--> Confirmed: {cfg['filename']} | Label: {cfg['label']} | Target Start: Year {cfg['start_year']}, Month {cfg['start_month']}")
    else:
        print(f"Warning: Explicitly specified file not found, skipping: {full_path}")

if not valid_files:
    raise FileNotFoundError("Error: None of the explicitly specified NetCDF files were found.")

# Collect distinct variables across all valid datasets to catch variables unique to Spin6
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

    # sharex=True and sharey=True enforce identical scale scopes globally across all panels
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 3.5 * nrows), dpi=150, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    legend_handles = []
    legend_labels = []

    for cfg in valid_files:
        with xr.open_dataset(cfg['path']) as ds:

            # Safely skip files (like Spin1 datasets) if they don't contain melt_rate_total
            if var_name not in ds.data_vars:
                continue

            raw_time = ds['Time'].values

            native_start_dec = raw_time[0].year + (raw_time[0].month - 1) / 12.0
            target_start_dec = cfg['start_year'] + (cfg['start_month'] - 1) / 12.0
            dynamic_offset = target_start_dec - native_start_dec

            time_coords = [
                (t.year + (t.month - 1) / 12.0) + dynamic_offset
                for t in raw_time
            ]

            data_slice = ds[var_name]

            for reg_idx, region in enumerate(regions):
                ax = axes_flat[reg_idx]

                if region in ds['region'].values:
                    reg_data = data_slice.sel(region=region)

                    # Extract raw numpy values for transformation
                    plot_values = reg_data.values

                    # Apply unit conversions depending on the targeted variable name
                    if 'melt_rate' in var_name:
                        # Conversion factor: kg/m^2/s -> m/year
                        plot_values = plot_values * sec_per_year / rho_fw
                    elif 'Uc' in var_name:
                        # Conversion factor: m/s -> cm/s
                        plot_values = plot_values * 100.0

                    line, = ax.plot(time_coords, plot_values,
                                    label=cfg['label'],
                                    linewidth=cfg.get('linewidth', 1.5),
                                    linestyle=cfg.get('linestyle', '-'),
                                    alpha=0.85,
                                    color=cfg.get('color', None))

                    # Track legend metrics only from the first subplot panel it encounters
                    if cfg.get('in_legend', True) and cfg['label'] not in legend_labels:
                        legend_handles.append(line)
                        legend_labels.append(cfg['label'])

    # Step 5: Formatting and beautifying the subplots
    for reg_idx, region in enumerate(regions):
        ax = axes_flat[reg_idx]
        ax.set_title(f"Region: {region}", fontsize=11, fontweight='normal', loc='left')
        ax.grid(True, linestyle='--', alpha=0.5)

        ax.ticklabel_format(style='plain', useOffset=False, axis='y')

        if reg_idx % ncols == 0:  # Only label leftmost y-axes
            if 'melt_rate' in var_name:
                ax.set_ylabel("Melt Rate Flux ($m/year$)", fontsize=10)
            elif 'Uc' in var_name:
                ax.set_ylabel("Friction Velocity ($cm/s$)", fontsize=10)
            elif 'Tc' in var_name:
                ax.set_ylabel("Thermal Driving ($\\Delta^\\circ C$)", fontsize=10)
            else:
                ax.set_ylabel("Value", fontsize=10)

        if reg_idx >= (num_regions - ncols):
            ax.set_xlabel("Adjusted Model Year", fontsize=10)

    # Clean up leftover empty subplot blocks
    for residual_idx in range(num_regions, len(axes_flat)):
        fig.delaxes(axes_flat[residual_idx])

    fig.suptitle(f"Regional Comparison Time Series: {var_name}", fontsize=16, fontweight='normal', y=0.99)

    if legend_handles:
        fig.legend(legend_handles, legend_labels, loc='upper center',
                   bbox_to_anchor=(0.5, 0.96), ncol=len(legend_handles), frameon=True, fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    save_path = os.path.join(output_dir, f"{var_name}.png")
    if opt_save == 1:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        fig.set_dpi(95)
        plt.show()
    plt.close(fig)

print(f"\nSuccess! All multi-panel comparison plots completed.")