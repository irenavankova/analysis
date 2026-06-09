#!/usr/bin/env python3

import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# --- User Configuration Block ---
# =============================================================================
# Specify a list of properties to loop through: 'temperature', 'salinity', 'Nsquare', 'rho'
TARGET_PROPERTIES = ['temperature', 'salinity', 'Nsquare', 'rho']

# Specify a list of time-series site configuration groups to process sequential loops for
TARGET_TS_NAMES = ['pts_berknerwest', 'pts_ronnedepr', 'pts_filchdepr', 'pts_ronnecenter', 'pts_shelfbreak', 'obs']

# Shading Toggle Switch: True will plot the time-series standard deviation bounds, False will hide it.
PLOT_STD_SHADING = True

# Order of resolutions to plot within each site subplot
resolutions = ['F8', 'F4', 'F2', 'F1']
res_colors = {'F1': 'black', 'F2': 'orange', 'F4': 'dodgerblue', 'F8': 'brown'}

# Define seasonal configurations: Name, month filter list, and label suffix
SEASON_CONFIGS = {
    'ANN': {'months': None, 'label': 'Annual'},
    'SO': {'months': [9, 10], 'label': 'Sep-Oct'},
    'DJF': {'months': [12, 1, 2], 'label': 'Dec-Feb'}
}

# --- Manual Limits for Profile Plots ---
vmin_temp, vmax_temp = -2.5, 0.0
vmin_salt, vmax_salt = 34.4, 35.0

# Metadata for Fields
param_meta = {
    'temperature': {
        'nc_var': 'timeMonthly_avg_activeTracers_temperature',
        'label': 'Temperature (°C)',
        'suffix': 'Temp',
        'xlim': None  # (vmin_temp, vmax_temp)
    },
    'salinity': {
        'nc_var': 'timeMonthly_avg_activeTracers_salinity',
        'label': 'Salinity (PSU)',
        'suffix': 'Sal',
        'xlim': None  # (vmin_salt, vmax_salt)
    },
    'Nsquare': {
        'nc_var': 'timeMonthly_avg_BruntVaisalaFreqTop',
        'label': 'Brunt-Väisälä Frequency ($s^{-2}$)',
        'suffix': 'Nsquare',
        'xlim': None
    },
    'rho': {
        'nc_var': 'timeMonthly_avg_potentialDensity',
        'label': r'Potential Density $\rho$ ($kg/m^3$)',
        'suffix': 'Rho',
        'xlim': None
    }
}


def get_field_data(ds_site, param, valid_layers_mask):
    if param in ds_site:
        dim_order = ('nVertLevels', 'Time') if 'nVertLevels' in ds_site[param].dims else (ds_site[param].dims[0],
                                                                                          'Time')
        return ds_site[param].transpose(*dim_order).values[valid_layers_mask, :]
    return None


# =============================================================================
# OUTER EXECUTION LOOPS: Iterating over configured keys
# =============================================================================
for ts_name in TARGET_TS_NAMES:

    input_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/{ts_name}'
    output_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_plots/vertical_profiles/{ts_name}'
    os.makedirs(output_dir, exist_ok=True)

    # Dynamic Reconstruction of Simulation File Target Structs
    simulation_cases = {
        'Spin6': {
            'F1': [f'{ts_name}_tseries_F1_Spin6p1.nc'],
            'F2': [f'{ts_name}_tseries_F2_Spin6p1.nc'],
            'F4': [f'{ts_name}_tseries_F4_Spin6p1.nc'],
            'F8': [f'{ts_name}_tseries_F8_Spin6p1.nc']
        }
    }

    all_possible_sites = []
    if ts_name == 'obs':
        all_possible_sites = [
            "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
            "R12", "R13", "R14", "F15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
            "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4", "FSW1"
        ]

    # --- Loading Block ---
    print("\n" + "=" * 85)
    print(f"LOADING SPIN6 DATASETS FOR GROUP: {ts_name}...")
    print("=" * 85)

    combined_datasets = {}

    for res in resolutions:
        spin6_files = simulation_cases['Spin6'].get(res, [])
        spin6_ds_list = []

        for idx, fname in enumerate(spin6_files):
            full_path = os.path.join(input_dir, fname)
            if os.path.exists(full_path):
                ds_part = xr.open_dataset(full_path)

                if not all_possible_sites and 'site' in ds_part:
                    all_possible_sites = list(ds_part['site'].values)
                    print(f" Found {len(all_possible_sites)} sites in: {fname}")

                spin6_ds_list.append(ds_part)

        if not spin6_ds_list:
            print(f"⚠️ Warning: Missing Spin6 files for resolution {res}. Skipping.")
            continue

        ds_spin6_full = xr.concat(spin6_ds_list, dim='Time')

        # Enforce continuous verification boundary checking
        if ds_spin6_full.sizes['Time'] < 48:
            print(
                f"❌ Error: Dataset for resolution {res} under group {ts_name} only has {ds_spin6_full.sizes['Time']} months.")
            print(f"   The requested analysis window (Years 2-4) requires at least 48 months of data.")
            sys.exit(1)

        ds_spin6_cut = ds_spin6_full.isel(Time=slice(12, 48))
        combined_datasets[res] = ds_spin6_cut

    # --- Seasonal Iteration Loop ---
    for season_id, config in SEASON_CONFIGS.items():
        print(f"\n--- Processing Season: {config['label']} ({season_id}) ---")

        # Apply seasonal slicing to combined datasets upfront
        seasonal_datasets = {}
        for res, ds in combined_datasets.items():
            if config['months'] is not None:
                # Filter indices where the Time dimension coordinate belongs to the target seasonal months
                seasonal_datasets[res] = ds.sel(Time=ds.Time.dt.month.isin(config['months']))
            else:
                seasonal_datasets[res] = ds

        # --- Plotting Block ---
        for target_prop in TARGET_PROPERTIES:
            meta = param_meta[target_prop]
            nc_var = meta['nc_var']

            print(f" Generating Unified Grid Plot Matrix -> Property: {meta['label']} ({season_id})")

            valid_sites = [s for s in all_possible_sites if any(
                res in seasonal_datasets and s in seasonal_datasets[res]['site'].values for res in resolutions)]
            num_sites = len(valid_sites)

            if num_sites == 0:
                print(
                    f"⚠️ Warning: No valid sites found under group {ts_name} for {target_prop} in season {season_id}. Skipping.")
                continue

            ncols = 4
            nrows = int(np.ceil(num_sites / ncols))

            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 5 * nrows),
                                     sharey=False)
            axes = axes.flatten()

            # Track global x bounds within this dynamic configuration loop
            x_min_fig = float('inf')
            x_max_fig = float('-inf')
            active_axes = []

            for s_idx, site_name in enumerate(valid_sites):
                ax = axes[s_idx]
                has_data = False

                z_min_site = float('inf')
                z_max_site = float('-inf')

                for res in resolutions:
                    if res not in seasonal_datasets:
                        continue

                    ds_res = seasonal_datasets[res]
                    if site_name not in ds_res['site'].values:
                        continue

                    ds_site = ds_res.sel(site=site_name)

                    # Build Vertical Grid Mesh using the first timestep of the filtered seasonal collection
                    ssh_t0 = float(
                        ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0
                    thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values
                    valid_layers_mask = ~np.isnan(thick_t0)
                    thick_valid = thick_t0[valid_layers_mask]

                    if len(thick_valid) == 0:
                        continue

                    z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)
                    z_centers = (z_edges[:-1] + z_edges[1:]) / 2.0

                    data_matrix = get_field_data(ds_site, nc_var, valid_layers_mask)
                    if data_matrix is None or data_matrix.size == 0:
                        continue

                    # Compute profile mean across Time dimension axis (restricted to chosen season)
                    mean_profile = np.nanmean(data_matrix, axis=1)

                    # Render mean profiles
                    ax.plot(mean_profile, z_centers, label=res, color=res_colors[res], linewidth=1.8, zorder=3)

                    # Check option switch to render standard deviation horizontal shading bounding region
                    if PLOT_STD_SHADING:
                        std_profile = np.nanstd(data_matrix, axis=1)
                        ax.fill_betweenx(z_centers, mean_profile - std_profile, mean_profile + std_profile,
                                         color=res_colors[res], alpha=0.15, zorder=2)

                    z_min_site = min(z_min_site, np.min(z_centers))
                    z_max_site = max(z_max_site, np.max(z_centers))
                    has_data = True

                if has_data:
                    ax.grid(True, linestyle=':', alpha=0.6)
                    ax.set_ylim(z_min_site, z_max_site)

                    if s_idx % ncols == 0:
                        ax.set_ylabel('Elevation z (m)', fontsize=11)
                    if s_idx >= (nrows - 1) * ncols:
                        ax.set_xlabel(meta['label'], fontsize=11)

                    if meta['xlim'] is not None:
                        x_min_fig, x_max_fig = meta['xlim']
                    else:
                        ax.relim()
                        ax.autoscale_view(scaley=False)
                        current_xmin, current_xmax = ax.get_xlim()
                        x_min_fig = min(x_min_fig, current_xmin)
                        x_max_fig = max(x_max_fig, current_xmax)

                    active_axes.append(ax)
                else:
                    ax.text(0.5, 0.5, "No data available", ha='center', va='center', transform=ax.transAxes,
                            color='grey')

            # --- Uniform X-Limit Synchronizer ---
            if active_axes and x_min_fig < x_max_fig:
                for ax in active_axes:
                    ax.set_xlim(x_min_fig, x_max_fig)

            # Hide excess empty subplot frames
            for empty_idx in range(num_sites, len(axes)):
                fig.delaxes(axes[empty_idx])

            axes[0].legend(title="Resolutions", loc='best')
            plt.tight_layout()

            # Append the season tracking identifier to the final filename structured block
            out_img_name = os.path.join(output_dir, f'Profiles_Spin6_Yrs2-4_{meta["suffix"]}_{season_id}.png')
            plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
            plt.close()
            print(f" -> Plot completely rendered and saved to: {out_img_name}")

    # Memory cleanup for loop dataset handles before advancing key groups
    for res, open_ds in combined_datasets.items():
        open_ds.close()

print("\nAll pipeline execution steps ran successfully through your target loops!")