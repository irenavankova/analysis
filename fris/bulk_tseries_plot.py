#!/usr/bin/env python3

import os
import glob
import xarray as xr
import matplotlib.pyplot as plt

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================
# Directory containing your reduced netCDF files ('.' means current directory)
data_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/'

# Define offsets for specific runs.
# If a file is found but its name isn't listed here, its offset defaults to 0.0.
year_offsets = {
    'F8': 0.0,
    'F4': 0.0,
    # Format: 'Unique_Filename_Part': offset_in_years
}

output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/bulk_tseries/compare_all/'
os.makedirs(output_dir, exist_ok=True)
# ==============================================================================

print(f"Scanning '{data_dir}' for time series netCDF files...")
file_pattern = os.path.join(data_dir, "bulk_tseries_*.nc")
found_files = sorted(glob.glob(file_pattern))

if not found_files:
    raise FileNotFoundError(f"Error: No files matching 'bulk_tseries_*.nc' found in '{data_dir}'")

# Automatically build the configuration list from discovered files
files_config = []
for fpath in found_files:
    fname = os.path.basename(fpath)

    # Generate a clean label by stripping 'bulk_tseries_' and '.nc'
    # e.g., 'bulk_tseries_F8_Spin6.nc' -> 'F8_Spin6'
    clean_label = fname.replace('bulk_tseries_', '').replace('.nc', '')

    # Check if this file or label matches any keys in our offset dictionary
    offset_value = 0.0
    for key, val in year_offsets.items():
        if key in fname:
            offset_value = val
            break

    files_config.append({
        'path': fpath,
        'label': clean_label,
        'offset': offset_value
    })
    print(f"--> Detected: {fname} | Assigned Label: {clean_label} | Offset: {offset_value} years")

# Step 1: Use the first file to determine the available variables and regions
with xr.open_dataset(files_config[0]['path']) as ds_ref:
    plot_vars = [v for v in ds_ref.data_vars if v not in ['Time', 'region']]
    regions = ds_ref['region'].values.tolist()

num_regions = len(regions)
print(f"\nProcessing {len(plot_vars)} variables across {num_regions} regions...")

# Determine the grid layout dynamically (2 columns wide)
ncols = 2
nrows = (num_regions + 1) // ncols

# Step 2: Loop over each variable (e.g., temp_max, salt_mean, rho_min)
for var_name in plot_vars:
    print(f"Generating multi-panel comparison plot for: {var_name}...")

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 3.5 * nrows), dpi=150, sharex=True)
    axes_flat = axes.flatten()

    legend_handles = []
    legend_labels = []

    # Step 3: Loop through each discovered file and plot its data onto the panels
    for cfg in files_config:
        with xr.open_dataset(cfg['path']) as ds:

            if var_name not in ds.data_vars:
                continue

            raw_time = ds['Time'].values

            # Apply decimal conversion + user-specified year offset
            time_coords = [
                (t.year + (t.month - 1) / 12.0) + cfg['offset']
                for t in raw_time
            ]

            data_slice = ds[var_name]

            # Step 4: Distribute data points to their respective regional subplots
            for reg_idx, region in enumerate(regions):
                ax = axes_flat[reg_idx]

                if region in ds['region'].values:
                    reg_data = data_slice.sel(region=region)

                    line, = ax.plot(time_coords, reg_data.values,
                                    label=cfg['label'], linewidth=1.5, alpha=0.85)

                    if reg_idx == 0:
                        legend_handles.append(line)
                        legend_labels.append(cfg['label'])

    # Step 5: Formatting and beautifying the subplots
    for reg_idx, region in enumerate(regions):
        ax = axes_flat[reg_idx]
        ax.set_title(f"Region: {region}", fontsize=11, fontweight='bold', loc='left')
        ax.grid(True, linestyle='--', alpha=0.5)

        if reg_idx % ncols == 0:  # Only label leftmost y-axes
            if 'temp' in var_name:
                ax.set_ylabel("Temperature (°C)", fontsize=10)
            elif 'salt' in var_name:
                ax.set_ylabel("Salinity (psu)", fontsize=10)
            elif 'rho' in var_name:
                ax.set_ylabel("Density (kg/m³)", fontsize=10)
            else:
                ax.set_ylabel("Value", fontsize=10)

        if reg_idx >= (num_regions - ncols):
            ax.set_xlabel("Adjusted Model Year (incl. Offset)", fontsize=10)

    # Clean up leftover empty subplot blocks
    for residual_idx in range(num_regions, len(axes_flat)):
        fig.delaxes(axes_flat[residual_idx])

    # Aesthetics and global legend positioning
    fig.suptitle(f"Regional Comparison Time Series: {var_name}", fontsize=16, fontweight='bold', y=0.99)
    fig.legend(legend_handles, legend_labels, loc='upper center',
               bbox_to_anchor=(0.5, 0.96), ncol=len(files_config), frameon=True, fontsize=11)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    save_path = os.path.join(output_dir, f"compare_{var_name}.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)

print(f"\nSuccess! All multi-panel comparison plots are saved in: '{output_dir}/'")