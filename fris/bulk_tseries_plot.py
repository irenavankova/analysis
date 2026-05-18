#!/usr/bin/env python3

import os
import glob
import xarray as xr
import matplotlib.pyplot as plt

Fnum = '8'
dx = f'F{Fnum}'
sec = 'Spin6'

# Define filename and output directory for plots
fpath = '/Users/ivankova/Desktop/Fris_hr/Fris_derived'
fsave = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/bulk_tseries'
nc_file = f'{fpath}/bulk_tseries_1by1_{dx}_{sec}.nc'
output_dir = f'{fsave}/plots_{dx}_{sec}'
os.makedirs(output_dir, exist_ok=True)

print(f"Opening data file: {nc_file}")
if not os.path.exists(nc_file):
    raise FileNotFoundError(f"Could not find the time series file: {nc_file}")

# Load the consolidated time series dataset
with xr.open_dataset(nc_file) as ds:
    # Identify variables to plot (exclude coordinate/dimension names)
    plot_vars = [v for v in ds.data_vars if v not in ['Time', 'region']]

    # Extract coordinates
    raw_time_coords = ds['Time'].values
    regions = ds['region'].values

    # --- FIX: Convert cftime to a decimal model year float array ---
    print("Converting climate calendar to continuous decimal model years...")
    time_coords = [
        t.year + (t.month - 1) / 12.0
        for t in raw_time_coords
    ]

    print(f"Found {len(plot_vars)} variables to plot over {len(regions)} regions.")

    # Loop over each reduced 3D variable (e.g., temp_max, salt_mean, rho_min)
    for var_name in plot_vars:
        print(f"Generating plot for: {var_name}...")

        # Initialize a crisp, clean figure
        fig, ax = plt.subplots(figsize=(11, 6), dpi=150)

        # Extract data slice for this variable: shape (Time, region)
        data_slice = ds[var_name]

        # Loop through regions and plot them as separate lines
        for region in regions:
            # Select data for the specific region
            reg_data = data_slice.sel(region=region)

            # Plot using the decimal model year axis
            ax.plot(time_coords, reg_data.values, label=str(region), linewidth=1.5)

        # Figure formatting and aesthetics
        ax.set_title(f"Time Series: {var_name} ({dx} - {sec})", fontsize=14, fontweight='normal', pad=15)
        ax.set_xlabel("Model Year", fontsize=12, labelpad=10)

        # Dynamically set y-axis labels based on the tracer type
        if 'temp' in var_name:
            ax.set_ylabel("Temperature (°C)", fontsize=12)
        elif 'salt' in var_name:
            ax.set_ylabel("Salinity (psu)", fontsize=12)
        elif 'rho' in var_name:
            ax.set_ylabel("Potential Density (kg/m³)", fontsize=12)
        else:
            ax.set_ylabel("Value", fontsize=12)

        ax.grid(True, linestyle='--', alpha=0.5)

        # Place legend nicely outside the plot frame
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), title="Regions", frameon=True)

        plt.tight_layout()

        # Define the output save path
        save_path = os.path.join(output_dir, f"{var_name}_{dx}_{sec}.png")
        plt.savefig(save_path, bbox_inches='tight')
        plt.close(fig)  # Free up memory immediately

print(f"Success! All plots have been saved into the directory: '{output_dir}/'")