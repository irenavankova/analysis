#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --- Configuration ---
input_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/obs_tseries'
output_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_plots/hovmoller'
os.makedirs(output_dir, exist_ok=True)

# Order of resolutions defining the vertical rows
resolutions = ['F1', 'F2', 'F4', 'F8']

# Chained file configurations by Simulation Case
simulation_cases = {
    'Spin1': {
        'F1': ['obs_tseries_F1_Spin1p1.nc', 'obs_tseries_F1_Spin1p2.nc', 'obs_tseries_F1_Spin1p3.nc'],
        'F2': ['obs_tseries_F2_Spin1p1.nc', 'obs_tseries_F2_Spin1p2.nc'],
        'F4': ['obs_tseries_F4_Spin1p1.nc'],
        'F8': ['obs_tseries_F8_Spin1p1.nc']
    },
    'Spin6': {
        'F1': ['obs_tseries_F1_Spin6p1.nc'],
        'F2': ['obs_tseries_F2_Spin6p1.nc'],
        'F4': ['obs_tseries_F4_Spin6p1.nc'],
        'F8': ['obs_tseries_F8_Spin6p1.nc']
    }
}

# Subplot columns meta configuration
parameters = ['temperature', 'salinity', 'flow_speed']
param_meta = {
    'temperature': {'label': 'Temperature (°C)', 'cmap': 'Spectral_r'},
    'salinity': {'label': 'Salinity (PSU)', 'cmap': 'viridis'},
    'flow_speed': {'label': 'Flow Speed (cm/s)', 'cmap': 'YlOrRd'}
}

# Max length to plot: 5 years * 12 months = 60 months
MAX_MONTHS = 60

# Global site coordinates reference list to loop through uniformly
all_possible_sites = [
    "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
    "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
    "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4",
    "FR5", "FR6", "FSW1"
]

# =============================================================================
# MAIN MULTI-DIMENSIONAL COORD LOOP
# =============================================================================
for sim_case, res_dict in simulation_cases.items():
    print("\n" + "=" * 85)
    print(f"LOADING DATASETS FOR SIMULATION CASE: {sim_case}")
    print("=" * 85)

    # Preload and merge resolution data for this case into a memory cache dict
    case_datasets = {}
    patch_times_dict = {}  # Store patch times in decimal years for each resolution

    for res in resolutions:
        file_list = res_dict.get(res, [])
        datasets_to_combine = []
        patch_months = []
        cumulative_months = 0

        for idx, fname in enumerate(file_list):
            full_path = os.path.join(input_dir, fname)
            if os.path.exists(full_path):
                ds_part = xr.open_dataset(full_path)
                datasets_to_combine.append(ds_part)

                # If this isn't the final segment, record where the next file patches in
                cumulative_months += ds_part.sizes['Time']
                if idx < len(file_list) - 1:
                    patch_months.append(cumulative_months)

        if not datasets_to_combine:
            print(f"⚠️ Warning: Missing files for {res}_{sim_case}. Row will be skipped.")
            continue

        # Merge time sequences and clip to max 5 years
        ds_full = xr.concat(datasets_to_combine, dim='Time')
        plot_months = min(ds_full.sizes['Time'], MAX_MONTHS)
        case_datasets[res] = ds_full.isel(Time=slice(0, plot_months))

        # Convert patch positions from month counts to decimal coordinate years
        res_patch_years = []
        for m in patch_months:
            if m < plot_months:  # Only display patches falling inside our 5-year window
                current_year = m // 12
                current_month_0indexed = m % 12
                decimal_year = current_year + (current_month_0indexed / 12.0)
                res_patch_years.append(decimal_year)

        patch_times_dict[res] = res_patch_years

    # Iterate site by site to construct the grid canvas matrix
    for site_name in all_possible_sites:

        # Lookahead validation check: does this site exist in at least one loaded resolution?
        site_exists_anywhere = False
        for res in resolutions:
            if res in case_datasets and site_name in case_datasets[res]['site'].values:
                site_exists_anywhere = True
                break

        if not site_exists_anywhere:
            continue

        print(f" -> Constructing Resolution vs Parameter Grid Canvas for Site: {site_name}...")

        # Initialize canvas grid figure: 4 rows (Resolutions) x 3 columns (Parameters)
        fig, axes = plt.subplots(nrows=len(resolutions), ncols=3, figsize=(22, 16), sharex='col')
        fig.suptitle(f'Resolution Comparison Hovmöller Profiles — Site: {site_name} ({sim_case})',
                     fontsize=18, fontweight='normal', y=0.96)

        # Fill the matrix grid row-by-row (Resolution), column-by-column (Parameter)
        for r_idx, res in enumerate(resolutions):

            # Label rows clearly on the left-most plots
            axes[r_idx, 0].annotate(f'{res}', xy=(0, 0.5), xytext=(-55, 0),
                                    xycoords='axes fraction', textcoords='offset points',
                                    size=14, weight='normal', ha='right', va='center', rotation=0)

            if res not in case_datasets:
                for c_idx in range(3):
                    axes[r_idx, c_idx].text(0.5, 0.5, f"No Data available for {res}", ha='center', va='center')
                continue

            ds_res = case_datasets[res]

            if site_name not in ds_res['site'].values:
                for c_idx in range(3):
                    axes[r_idx, c_idx].text(0.5, 0.5, f"Site missing in {res}", ha='center', va='center')
                continue

            ds_site = ds_res.sel(site=site_name)
            plot_months = ds_site.sizes['Time']

            # --- Coordinate Configurations (Time / X-Axis) ---
            time_coords = []
            tick_positions = []
            tick_labels = []

            for i in range(plot_months):
                current_year = i // 12
                current_month_0indexed = i % 12
                decimal_year = current_year + (current_month_0indexed / 12.0)
                time_coords.append(decimal_year)

                # Only store ticks and labels for Year X, Month 1 (index 0)
                if current_month_0indexed == 0:
                    tick_positions.append(decimal_year)
                    tick_labels.append(f"Yr {current_year} M1")

            time_coords = np.array(time_coords)

            if plot_months > 1:
                dx = time_coords[1] - time_coords[0]
            else:
                dx = 1.0 / 12.0
            x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

            # --- Coordinate Configurations (Depth / Y-Axis) ---
            ssh_t0 = float(
                ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0
            thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values
            valid_layers_mask = ~np.isnan(thick_t0)
            thick_valid = thick_t0[valid_layers_mask]

            if len(thick_valid) == 0:
                for c_idx in range(3):
                    axes[r_idx, c_idx].text(0.5, 0.5, "No active ocean layers", ha='center', va='center')
                continue

            z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)

            # --- Field Derivation & Grid Plotting ---
            for c_idx, param in enumerate(parameters):
                ax = axes[r_idx, c_idx]
                meta = param_meta[param]

                # Dynamic variable isolation
                if param == 'temperature':
                    vname = 'timeMonthly_avg_activeTracers_temperature'
                    data_matrix = ds_site[vname].transpose('nVertLevels', 'Time').values[
                        valid_layers_mask, :] if vname in ds_site else None
                elif param == 'salinity':
                    vname = 'timeMonthly_avg_activeTracers_salinity'
                    data_matrix = ds_site[vname].transpose('nVertLevels', 'Time').values[
                        valid_layers_mask, :] if vname in ds_site else None
                elif param == 'flow_speed':
                    v_mer = 'timeMonthly_avg_velocityMeridional'
                    v_zon = 'timeMonthly_avg_velocityZonal'
                    if v_mer in ds_site and v_zon in ds_site:
                        m_matrix = ds_site[v_mer].transpose('nVertLevels', 'Time').values[valid_layers_mask, :]
                        z_matrix = ds_site[v_zon].transpose('nVertLevels', 'Time').values[valid_layers_mask, :]
                        # FIX: Apply scaling multiplier outside the square root operation to preserve velocity magnitudes
                        data_matrix = np.sqrt(m_matrix ** 2 + z_matrix ** 2) * 100
                    else:
                        data_matrix = None

                if data_matrix is None:
                    ax.text(0.5, 0.5, "Missing Fields", ha='center', va='center', color='red')
                    continue

                # Render pcolormesh field
                mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'], edgecolors='none')
                ax.grid(True, linestyle=':', alpha=0.3, color='black')

                # --- Draw Vertical Patch Boundary Lines ---
                # Retrieve recorded timestamps for this resolution row and plot them across columns
                for patch_time in patch_times_dict.get(res, []):
                    ax.axvline(x=patch_time, color='black', linestyle='-', linewidth=1.5, zorder=3)

                # Formatting Labels for Top and Bottom Bounds
                if r_idx == 0:
                    ax.set_title(meta['label'], fontsize=14, pad=12, fontweight='normal')
                if c_idx == 0:
                    ax.set_ylabel('Elevation z (m)', fontsize=11)

                # Set up shared axis properties for bottom row elements
                if r_idx == len(resolutions) - 1:
                    ax.set_xticks(tick_positions)
                    ax.set_xticklabels(tick_labels, rotation=0, ha='center', fontsize=10)
                    ax.set_xlabel('Model Year', fontsize=12, labelpad=10)

                ax.set_xlim(x_edges[0], x_edges[-1])

                # Add standalone colorbar tightly paired against each panel plot
                cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, aspect=18)
                cbar.ax.tick_params(labelsize=9)

        # Space out subplots cleanly and export
        plt.subplots_adjust(hspace=0.20, wspace=0.18, bottom=0.08, top=0.89, left=0.12, right=0.95)
        out_img_name = os.path.join(output_dir, f'hovmoller_grid_{site_name}_{sim_case}_5yr.png')
        plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
        plt.close()

    # Clean cache objects out of memory loop
    for res, open_ds in case_datasets.items():
        open_ds.close()

print("\n Combined resolution/parameter grid comparison plots completed successfully!")