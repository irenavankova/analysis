#!/usr/bin/env python3

import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --- Configuration ---
ts_name = 'pts_berknerwest'
ts_name = 'pts_ronnedepr'
ts_name = 'pts_filchdepr'
ts_name = 'pts_ronnecenter'
ts_name = 'pts_shelfbreak'
ts_name = 'obs'





input_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/{ts_name}'
output_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_plots/hovmoller/TS/{ts_name}'
os.makedirs(output_dir, exist_ok=True)

# Order of resolutions defining the vertical rows
resolutions = ['F1', 'F2', 'F4', 'F8']

# --- Manual Colorbar Limits (Edit Here) ---
vmin_temp, vmax_temp = -2.5, -1.3  # Min/Max for Temperature (°C)
vmin_salt, vmax_salt = 34.4, 35.0  # Min/Max for Salinity (PSU)

# Chained file configurations by Simulation Case
simulation_cases = {
    'Spin1': {
        'F1': [f'{ts_name}_tseries_F1_Spin1p1.nc', f'{ts_name}_tseries_F1_Spin1p2.nc',
               f'{ts_name}_tseries_F1_Spin1p3.nc'],
        'F2': [f'{ts_name}_tseries_F2_Spin1p1.nc', f'{ts_name}_tseries_F2_Spin1p2.nc'],
        'F4': [f'{ts_name}_tseries_F4_Spin1p1.nc'],
        'F8': [f'{ts_name}_tseries_F8_Spin1p1.nc']
    },
    'Spin6': {
        'F1': [f'{ts_name}_tseries_F1_Spin6p1.nc'],
        'F2': [f'{ts_name}_tseries_F2_Spin6p1.nc'],
        'F4': [f'{ts_name}_tseries_F4_Spin6p1.nc'],
        'F8': [f'{ts_name}_tseries_F8_Spin6p1.nc']
    }
}

# Metadata for Targeted Fields
param_meta = {
    'timeMonthly_avg_activeTracers_temperature': {
        'label': 'Temperature (°C)',
        'cmap': 'Spectral_r',
        'vmin': vmin_temp,
        'vmax': vmax_temp,
        'suffix': 'Temp',
        'contour_val': -1.9
    },
    'timeMonthly_avg_activeTracers_salinity': {
        'label': 'Salinity (PSU)',
        'cmap': 'viridis',
        'vmin': vmin_salt,
        'vmax': vmax_salt,
        'suffix': 'Sal',
        'contour_val': 34.9
    }
}

MAX_MONTHS = 60  # 5 years for each simulation case

all_possible_sites = []
if ts_name == 'obs':
    all_possible_sites = [
        "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
        "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
        "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4", "FSW1"
    ]


def get_field_data(ds_site, param, valid_layers_mask):
    if param in ds_site:
        dim_order = ('nVertLevels', 'Time') if 'nVertLevels' in ds_site[param].dims else (ds_site[param].dims[0],
                                                                                          'Time')
        return ds_site[param].transpose(*dim_order).values[valid_layers_mask, :]
    return None


# =============================================================================
# DATA LOADING BLOCK: Combine Spin1 & Spin6 Sequentially
# =============================================================================
print("\n" + "=" * 85)
print("LOADING AND COMBINING SPIN1 & SPIN6 DATASETS...")
print("=" * 85)

combined_datasets = {}
patch_times_dict = {}

for res in resolutions:
    # --- Load Spin1 (First 5 Years) ---
    spin1_files = simulation_cases['Spin1'].get(res, [])
    spin1_ds_list = []
    spin1_months_cumulative = 0
    res_patch_years = []

    for idx, fname in enumerate(spin1_files):
        full_path = os.path.join(input_dir, fname)
        if os.path.exists(full_path):
            ds_part = xr.open_dataset(full_path)

            if not all_possible_sites and 'site' in ds_part:
                all_possible_sites = list(ds_part['site'].values)
                print(f" Found {len(all_possible_sites)} sites in: {fname}")

            spin1_ds_list.append(ds_part)
            spin1_months_cumulative += ds_part.sizes['Time']
            if idx < len(spin1_files) - 1 and spin1_months_cumulative < MAX_MONTHS:
                res_patch_years.append(spin1_months_cumulative / 12.0)

    if not spin1_ds_list:
        print(f"⚠️ Warning: Missing Spin1 files for {res}. Skipping.")
        continue

    ds_spin1_full = xr.concat(spin1_ds_list, dim='Time')
    ds_spin1_cut = ds_spin1_full.isel(Time=slice(0, min(ds_spin1_full.sizes['Time'], MAX_MONTHS)))
    actual_spin1_months = ds_spin1_cut.sizes['Time']

    # --- Load Spin6 (First 5 Years) ---
    spin6_files = simulation_cases['Spin6'].get(res, [])
    spin6_ds_list = []
    spin6_months_cumulative = 0

    for idx, fname in enumerate(spin6_files):
        full_path = os.path.join(input_dir, fname)
        if os.path.exists(full_path):
            ds_part = xr.open_dataset(full_path)
            spin6_ds_list.append(ds_part)

            spin6_months_cumulative += ds_part.sizes['Time']
            if idx < len(spin6_files) - 1 and spin6_months_cumulative < MAX_MONTHS:
                res_patch_years.append((actual_spin1_months + spin6_months_cumulative) / 12.0)

    if not spin6_ds_list:
        print(f"⚠️ Warning: Missing Spin6 files for {res}. Plot will only show Spin1.")
        combined_datasets[res] = ds_spin1_cut
        patch_times_dict[res] = res_patch_years
        continue

    ds_spin6_full = xr.concat(spin6_ds_list, dim='Time')
    ds_spin6_cut = ds_spin6_full.isel(Time=slice(0, min(ds_spin6_full.sizes['Time'], MAX_MONTHS)))

    # Concatenate along continuous timeline coordinate axis
    ds_combined = xr.concat([ds_spin1_cut, ds_spin6_cut], dim='Time', data_vars='minimal', coords='minimal',
                            compat='override')
    combined_datasets[res] = ds_combined
    patch_times_dict[res] = res_patch_years

# =============================================================================
# MAIN PLOTTING LOOP
# =============================================================================
for site_name in all_possible_sites:

    # Verify site exists within the processed timelines
    site_exists_anywhere = any(
        res in combined_datasets and site_name in combined_datasets[res]['site'].values for res in resolutions)
    if not site_exists_anywhere:
        continue

    for param, meta in param_meta.items():
        print(f" -> Site: {site_name} | Processing Profile: {meta['label']}...")

        # Setup 4 rows by 1 column figure with shared x-axis mapping
        fig, axes = plt.subplots(nrows=len(resolutions), ncols=1, figsize=(8, 16), sharex=True)

        fig.suptitle(f'Resolution Comparison Hovmöller Profiles — Site: {site_name} ({meta["label"]})',
                     fontsize=15, fontweight='normal', y=0.94)

        for r_idx, res in enumerate(resolutions):
            ax = axes[r_idx]

            # Label vertical rows with resolution values on the left
            ax.annotate(f'{res}', xy=(0, 0.5), xytext=(-55, 0),
                        xycoords='axes fraction', textcoords='offset points',
                        size=14, weight='normal', ha='right', va='center')

            if res not in combined_datasets:
                ax.text(0.5, 0.5, f"No Data available for {res}", ha='center', va='center')
                continue

            ds_res = combined_datasets[res]
            if site_name not in ds_res['site'].values:
                ax.text(0.5, 0.5, f"Site missing in {res}", ha='center', va='center')
                continue

            ds_site = ds_res.sel(site=site_name)
            total_months = ds_site.sizes['Time']

            # Build Continuous Combined Time Coordinates Map
            time_coords, tick_positions, tick_labels = [], [], []
            for i in range(total_months):
                decimal_year = i / 12.0
                time_coords.append(decimal_year)

                if i % 12 == 0:
                    tick_positions.append(decimal_year)
                    if i < MAX_MONTHS:
                        tick_labels.append(f"S1 Yr {i // 12 + 1}")
                    else:
                        tick_labels.append(f"S6 Yr {(i - MAX_MONTHS) // 12 + 1}")

            time_coords = np.array(time_coords)
            dx = 1.0 / 12.0
            x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

            # Build Vertical Grid Mesh
            ssh_t0 = float(
                ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0
            thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values
            valid_layers_mask = ~np.isnan(thick_t0)
            thick_valid = thick_t0[valid_layers_mask]

            if len(thick_valid) == 0:
                ax.text(0.5, 0.5, "No active ocean layers", ha='center', va='center')
                continue

            z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)

            # Extract 2D grid matrix and mesh plot using manual min/max bounds
            data_matrix = get_field_data(ds_site, param, valid_layers_mask)
            if data_matrix is None:
                ax.text(0.5, 0.5, "Missing Fields", ha='center', va='center', color='red')
                continue

            mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'],
                                 vmin=meta['vmin'], vmax=meta['vmax'], edgecolors='none')

            # --- ADDED: Black Contour Lines for Specific Thresholds ---
            # Calculated using cell midpoints for precise tracking mapping
            z_centers = (z_edges[:-1] + z_edges[1:]) / 2.0
            if meta['contour_val'] >= np.nanmin(data_matrix) and meta['contour_val'] <= np.nanmax(data_matrix):
                contours = ax.contour(time_coords, z_centers, data_matrix,
                                      levels=[meta['contour_val']], colors='black',
                                      linewidths=1.2, linestyles='-')
                # Optional: Uncomment if you want text labels on the contour lines
                # ax.clabel(contours, inline=True, fmt=f"{meta['contour_val']}", fontsize=8)

            ax.grid(True, linestyle=':', alpha=0.3, color='black')

            # Magenta marker dividing Spin1 and Spin6 at 5 years
            if total_months > MAX_MONTHS:
                ax.axvline(x=MAX_MONTHS / 12.0, color='magenta', linestyle='-.', linewidth=1.0, zorder=4)

            # Segment breakdown patch dividers
            for patch_time in patch_times_dict.get(res, []):
                ax.axvline(x=patch_time, color='black', linestyle='--', linewidth=1.0, zorder=3)

            ax.set_ylabel('Elevation z (m)', fontsize=11)

            # Apply x-axis timeline designations to the bottom subplot panel
            if r_idx == len(resolutions) - 1:
                ax.set_xticks(tick_positions)
                ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=9)
                ax.set_xlabel('Simulation Run Timeline', fontsize=12, labelpad=10)

            ax.set_xlim(x_edges[0], x_edges[-1])
            cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, aspect=18)
            cbar.ax.tick_params(labelsize=9)

        plt.subplots_adjust(hspace=0.25, bottom=0.12, top=0.88, left=0.15, right=0.92)

        # Dynamically find and create exact subfolder paths safely
        out_img_name = os.path.join(output_dir, meta["suffix"], f'{site_name}_{meta["suffix"]}.png')
        os.makedirs(os.path.dirname(out_img_name), exist_ok=True)

        plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
        plt.close()

for res, open_ds in combined_datasets.items():
    open_ds.close()

print("\n Combined Spin1 & Spin6 Hovmöller profile processing completed successfully!")