#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --- Configuration ---
input_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/obs_tseries'
output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/hovmoller/hov_singleF'
os.makedirs(output_dir, exist_ok=True)

# Define your target groupings to ensure precise sequential chaining
simulations_groups = {
    'F1_Spin1': ['obs_tseries_F1_Spin1p1.nc', 'obs_tseries_F1_Spin1p2.nc', 'obs_tseries_F1_Spin1p3.nc'],
    'F1_Spin6': ['obs_tseries_F1_Spin6p1.nc'],
    'F2_Spin1': ['obs_tseries_F2_Spin1p1.nc', 'obs_tseries_F2_Spin1p2.nc'],
    'F2_Spin6': ['obs_tseries_F2_Spin6p1.nc'],
    'F4_Spin1': ['obs_tseries_F4_Spin1p1.nc'],
    'F4_Spin6': ['obs_tseries_F4_Spin6p1.nc'],
    'F8_Spin1': ['obs_tseries_F8_Spin1p1.nc'],
    'F8_Spin6': ['obs_tseries_F8_Spin6p1.nc']
}

fields_to_plot = [
    'timeMonthly_avg_velocityMeridional',
    'timeMonthly_avg_velocityZonal',
    'timeMonthly_avg_activeTracers_temperature',
    'timeMonthly_avg_activeTracers_salinity',
    'timeMonthly_avg_potentialDensity'
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

# Max length to plot: 5 years * 12 months = 60 months
MAX_MONTHS = 60

# =============================================================================
# MAIN LOGIC BLOCK: CHAINING, SPLICING & PLOTTING
# =============================================================================
for group_id, file_list in simulations_groups.items():
    f_res, sim_case = group_id.split('_')  # e.g., 'F1', 'Spin1'

    # 1. Filter out missing files and read available segments
    datasets_to_combine = []
    for fname in file_list:
        full_path = os.path.join(input_dir, fname)
        if os.path.exists(full_path):
            datasets_to_combine.append(xr.open_dataset(full_path))
        else:
            print(f"⚠️ Warning: Component file '{fname}' not found. Skipping segment.")

    if not datasets_to_combine:
        print(f"❌ Skipping configuration '{group_id}': No data files found.")
        continue

    print("\n" + "=" * 75)
    print(f"PROCESSING GROUP: {f_res} | {sim_case} (Chaining {len(datasets_to_combine)} parts)")
    print("=" * 75)

    # 2. Concatenate part sequences along the native 'Time' dimension
    ds_full = xr.concat(datasets_to_combine, dim='Time')

    # 3. Truncate dynamically to a maximum of 5 years (60 months)
    actual_months = ds_full.sizes['Time']
    #plot_months = min(actual_months, MAX_MONTHS)
    plot_months = actual_months

    # Slice the combined dataset to the target window
    ds = ds_full.isel(Time=slice(0, plot_months))
    print(f" -> Combined time series contains {actual_months} months. Plotting first {plot_months} months.")

    # 4. Generate Adjusted Coordinate Time System (relative decimal years from 0)
    time_coords = []
    time_labels = []
    for i in range(plot_months):
        current_year = i // 12
        current_month_0indexed = i % 12
        decimal_year = current_year + (current_month_0indexed / 12.0)

        time_coords.append(decimal_year)
        time_labels.append(f"Yr {current_year} M{current_month_0indexed + 1}")

    time_coords = np.array(time_coords)

    # Build 1D continuous boundary pixel edges for pcolormesh along Time
    if plot_months > 1:
        dx = time_coords[1] - time_coords[0]
    else:
        dx = 1.0 / 12.0
    x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

    sites = ds['site'].values

    # 5. Iterative Generation Loop Over Individual Grid Sites
    for site_name in sites:
        print(f"    -> Plotting site: {site_name}...")
        ds_site = ds.sel(site=site_name)

        # --- Reconstruct Physical Vertical Elevation (z) ---
        ssh_t0 = float(ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0
        thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values

        # Mask out native empty NaN layers safely
        valid_layers_mask = ~np.isnan(thick_t0)
        thick_valid = thick_t0[valid_layers_mask]

        if len(thick_valid) == 0:
            print(f"       ⚠️ Skipping: No valid active water column layers at {site_name}.")
            continue

        # Build cellular continuous boundary interface edges downward from SSH elevation
        z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)

        # Initialize subplot grid
        fig, axes = plt.subplots(nrows=len(fields_to_plot), ncols=1, figsize=(13, 18), sharex=True)
        fig.suptitle(f'Hovmöller Profiles — Site: {site_name} ({f_res} {sim_case} - First 5 Years Max)',
                     fontsize=15, fontweight='bold', y=0.95)

        for idx, field in enumerate(fields_to_plot):
            ax = axes[idx]

            if field not in ds_site:
                ax.text(0.5, 0.5, f"Variable '{field}' missing from dataset",
                        ha='center', va='center', transform=ax.transAxes, color='red')
                continue

            # Extract 2D matrix shape: (Time, nVertLevels) -> transpose to (nVertLevels, Time)
            full_matrix = ds_site[field].transpose('nVertLevels', 'Time').values

            # Filter active rows matching non-NaN cell coordinates
            data_matrix = full_matrix[valid_layers_mask, :]

            meta = field_meta.get(field, {'label': field, 'cmap': 'viridis'})

            # Scale divergent velocity metrics symmetrically around zero
            if 'velocity' in field:
                vmax = max(abs(np.nanmin(data_matrix)), abs(np.nanmax(data_matrix)))
                vmin = -vmax if vmax != 0 else -1
                vmax = vmax if vmax != 0 else 1
            else:
                vmin, vmax = None, None

            # Render data matrix explicitly over coordinates
            mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'], vmin=vmin, vmax=vmax,
                                 edgecolors='none')

            ax.set_ylabel('Elevation z (m)', fontsize=11)
            ax.grid(True, linestyle=':', alpha=0.3, color='black')

            # Build standardized uniform colorbars
            cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.015, aspect=12)
            cbar.set_label(meta['label'], fontsize=10, weight='bold')

        # Format and label the shared bottom X-axis
        ax_bottom = axes[-1]
        step = max(1, plot_months // 12)  # Thin out tick marks dynamically based on timeframe

        ax_bottom.set_xticks(time_coords[::step])
        ax_bottom.set_xticklabels([time_labels[i] for i in range(0, plot_months, step)], rotation=45, ha='right',
                                  fontsize=9)
        ax_bottom.set_xlabel('Adjusted Model Year', fontsize=12, labelpad=10)
        ax_bottom.set_xlim(x_edges[0], x_edges[-1])

        plt.subplots_adjust(hspace=0.28, bottom=0.10, top=0.91, left=0.10, right=0.92)

        # Clean file exports using the unified name identifiers
        out_img_name = os.path.join(output_dir, f'hovmoller_{site_name}_{f_res}_{sim_case}_5yr.png')
        plt.savefig(out_img_name, dpi=250, bbox_inches='tight')
        plt.close()

    # Close component references to manage memory
    for open_ds in datasets_to_combine:
        open_ds.close()

print("\n Analysis, file-group chaining, and visualization complete!")