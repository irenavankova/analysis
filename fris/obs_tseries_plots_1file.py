#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --- Configuration ---
input_file = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/obs_tseries/obs_tseries_F8_Spin1p1.nc'
output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/hovmoller'
os.makedirs(output_dir, exist_ok=True)

# Monthly offsets matching parameters for Spin1p1
start_year = 0
start_month = 1

fields_to_plot = [
    'timeMonthly_avg_velocityMeridional',
    'timeMonthly_avg_velocityZonal',
    'timeMonthly_avg_activeTracers_temperature',
    'timeMonthly_avg_activeTracers_salinity',
    'timeMonthly_avg_potentialDensity',
    'timeMonthly_avg_layerThickness'
]

# Clean display names and colormaps for oceanographic variables
field_meta = {
    'timeMonthly_avg_velocityMeridional': {'label': 'Meridional Velocity (m/s)',
                                           'cmap': 'cmr.ice_r' if 'cmr' in globals() else 'RdBu_r'},
    'timeMonthly_avg_velocityZonal': {'label': 'Zonal Velocity (m/s)',
                                      'cmap': 'cmr.ice_r' if 'cmr' in globals() else 'RdBu_r'},
    'timeMonthly_avg_activeTracers_temperature': {'label': 'Temperature (°C)', 'cmap': 'Spectral_r'},
    'timeMonthly_avg_activeTracers_salinity': {'label': 'Salinity (psu)', 'cmap': 'viridis'},
    'timeMonthly_avg_potentialDensity': {'label': 'Potential Density (kg/m³)', 'cmap': 'cividis'},
    'timeMonthly_avg_layerThickness': {'label': 'Layer Thickness (m)', 'cmap': 'GnBu'}
}

# --- Load Dataset ---
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Target profile dataset '{input_file}' not found. Please verify the path.")

print(f"Opening {input_file}...")
ds = xr.open_dataset(input_file)

# --- Calculate Adjusted Time Coordinates ---
num_timesteps = ds.sizes['Time']
time_coords = []
time_labels = []

# Generate relative decimal coordinates and string values starting from Year 0
for i in range(num_timesteps):
    total_months = (start_month - 1) + i
    current_year = start_year + (total_months // 12)
    current_month_0indexed = total_months % 12
    decimal_year = current_year + (current_month_0indexed / 12.0)

    time_coords.append(decimal_year)
    time_labels.append(f"Yr {current_year} M{current_month_0indexed + 1}")

time_coords = np.array(time_coords)

# Construct 1D boundary edges for the X-axis grid in pcolormesh
if num_timesteps > 1:
    dx = time_coords[1] - time_coords[0]
else:
    dx = 1.0 / 12.0
x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

sites = ds['site'].values
print(f"Found {len(sites)} sites to plot: {list(sites)}")

# --- Generation Loop ---
for site_name in sites:
    print(f"Generating Hovmöller multi-plot for Site: {site_name}...")

    # Isolate site profile data
    ds_site = ds.sel(site=site_name)

    # --- Vertical Coordinate (z) Reconstruction ---
    # Extract Sea Surface Height (ssh) at Time=0, fallback to 0.0 if missing
    ssh_t0 = float(ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0

    # Extract layer thicknesses from the first time step
    thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values

    # Mask out native NaN layers (bedrock / ice boundaries)
    valid_layers_mask = ~np.isnan(thick_t0)
    thick_valid = thick_t0[valid_layers_mask]

    if len(thick_valid) == 0:
        print(f"⚠️ Warning: No valid active water column layers for Site {site_name}. Skipping.")
        continue

    # Build continuous cell-interface z-edges downwards from SSH
    # z_edges start at ssh and progress more negative into depth: e.g. [0.5, -9.5, -21.0...]
    z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)

    # Initialize a stacked figure layout (6 subplots, sharing the same Time X-axis)
    fig, axes = plt.subplots(nrows=len(fields_to_plot), ncols=1, figsize=(13, 18), sharex=True)
    fig.suptitle(f'Hovmöller Profiles — Site: {site_name}', fontsize=16, fontweight='bold', y=0.95)

    for idx, field in enumerate(fields_to_plot):
        ax = axes[idx]

        if field not in ds_site:
            ax.text(0.5, 0.5, f"Variable '{field}' missing from dataset",
                    ha='center', va='center', transform=ax.transAxes, color='red')
            continue

        # Extract 2D grid matrix shape: (Time, nVertLevels) -> transpose to (nVertLevels, Time)
        full_matrix = ds_site[field].transpose('nVertLevels', 'Time').values

        # Filter rows corresponding to active valid vertical coordinates only
        data_matrix = full_matrix[valid_layers_mask, :]

        # Pull specific label details
        meta = field_meta.get(field, {'label': field, 'cmap': 'viridis'})

        # Center divergent colormaps symmetrically around 0 for velocities
        if 'velocity' in field:
            vmax = max(abs(np.nanmin(data_matrix)), abs(np.nanmax(data_matrix)))
            vmin = -vmax if vmax != 0 else -1
            vmax = vmax if vmax != 0 else 1
        else:
            vmin, vmax = None, None

        # Render grid directly onto physical vertical z-coordinates
        mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'], vmin=vmin, vmax=vmax, edgecolors='none')

        # Styling adjustments
        ax.set_ylabel('Elevation z (m)', fontsize=11)
        ax.grid(True, linestyle=':', alpha=0.3, color='black')

        # Add synchronized colorbars alongside each field plot
        cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.015, aspect=12)
        cbar.set_label(meta['label'], fontsize=10, weight='bold')

    # Configure Shared X-Axis labels using manual coordinate values
    ax_bottom = axes[-1]

    # Calculate sensible tick intervals depending on how long the time series is
    step = max(1, num_timesteps // 12)

    ax_bottom.set_xticks(time_coords[::step])
    ax_bottom.set_xticklabels([time_labels[i] for i in range(0, num_timesteps, step)], rotation=45, ha='right',
                              fontsize=9)
    ax_bottom.set_xlabel('Adjusted Model Year', fontsize=12, labelpad=10)
    ax_bottom.set_xlim(x_edges[0], x_edges[-1])

    # Clean spacing layouts
    plt.subplots_adjust(hspace=0.28, bottom=0.10, top=0.91, left=0.10, right=0.92)

    # Save chart out to high-res file inside the specified output path
    out_img_name = os.path.join(output_dir, f'hovmoller_{site_name}_F8_Spin1p1.png')
    plt.savefig(out_img_name, dpi=250, bbox_inches='tight')
    plt.close()

    print(f"Successfully exported: {out_img_name}")

print("\nAll Hovmöller plotting streams completed successfully.")