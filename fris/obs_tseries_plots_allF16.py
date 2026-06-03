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

# --- Dynamic Plotting Groups Mapping ---
plotting_groups = {
    'TSU': ['timeMonthly_avg_activeTracers_temperature', 'timeMonthly_avg_activeTracers_salinity', 'flow_speed'],
    'Density': ['timeMonthly_avg_potentialDensity', 'timeMonthly_avg_BruntVaisalaFreqTop',
                'timeMonthly_avg_layerThickness'],
    'Tadv': ['timeMonthly_avg_activeTracersTend_temperatureTend',
             'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency',
             'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency'],
    'Tmix': ['timeMonthly_avg_activeTracerVertMixTendency_temperatureVertMixTendency',
             'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency',
             'timeMonthly_avg_activeTracerNonLocalTendency_temperatureNonLocalTendency'],
    'Sadv': ['timeMonthly_avg_activeTracersTend_salinityTend',
             'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_salinityHorizontalAdvectionTendency',
             'timeMonthly_avg_activeTracerVerticalAdvectionTendency_salinityVerticalAdvectionTendency'],
    'Smix': ['timeMonthly_avg_activeTracerVertMixTendency_salinityVertMixTendency',
             'timeMonthly_avg_activeTracerHorMixTendency_salinityHorMixTendency',
             'timeMonthly_avg_activeTracerNonLocalTendency_salinityNonLocalTendency']
}

# Metadata Dictionary mapping netCDF variable names to clear titles & logical colormaps
param_meta = {
    'timeMonthly_avg_activeTracers_temperature': {'label': 'Temperature (°C)', 'cmap': 'Spectral_r'},
    'timeMonthly_avg_activeTracers_salinity': {'label': 'Salinity (PSU)', 'cmap': 'viridis'},
    'flow_speed': {'label': 'Flow Speed (cm/s)', 'cmap': 'YlOrRd'},
    'timeMonthly_avg_potentialDensity': {'label': 'Potential Density (kg/m³)', 'cmap': 'cividis'},
    'timeMonthly_avg_BruntVaisalaFreqTop': {'label': 'Brunt-Väisälä Freq (s⁻²)', 'cmap': 'magma'},
    'timeMonthly_avg_layerThickness': {'label': 'Layer Thickness (m)', 'cmap': 'Blues'},

    # Temperature Tendencies
    'timeMonthly_avg_activeTracersTend_temperatureTend': {'label': 'Total Temp Tendency (°C/s)', 'cmap': 'RdBu_r'},
    'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency': {
        'label': 'Temp Horiz Advection Tend (°C/s)', 'cmap': 'RdBu_r'},
    'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency': {
        'label': 'Temp Vert Advection Tend (°C/s)', 'cmap': 'RdBu_r'},
    'timeMonthly_avg_activeTracerVertMixTendency_temperatureVertMixTendency': {'label': 'Temp Vert Mix Tend (°C/s)',
                                                                               'cmap': 'RdBu_r'},
    'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency': {'label': 'Temp Horiz Mix Tend (°C/s)',
                                                                             'cmap': 'RdBu_r'},
    'timeMonthly_avg_activeTracerNonLocalTendency_temperatureNonLocalTendency': {'label': 'Temp Non-local Tend (°C/s)',
                                                                                 'cmap': 'RdBu_r'},

    # Salinity Tendencies
    'timeMonthly_avg_activeTracersTend_salinityTend': {'label': 'Total Salt Tendency (PSU/s)', 'cmap': 'BrBG_r'},
    'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_salinityHorizontalAdvectionTendency': {
        'label': 'Salt Horiz Advection Tend (PSU/s)', 'cmap': 'BrBG_r'},
    'timeMonthly_avg_activeTracerVerticalAdvectionTendency_salinityVerticalAdvectionTendency': {
        'label': 'Salt Vert Advection Tend (PSU/s)', 'cmap': 'BrBG_r'},
    'timeMonthly_avg_activeTracerVertMixTendency_salinityVertMixTendency': {'label': 'Salt Vert Mix Tend (PSU/s)',
                                                                            'cmap': 'BrBG_r'},
    'timeMonthly_avg_activeTracerHorMixTendency_salinityHorMixTendency': {'label': 'Salt Horiz Mix Tend (PSU/s)',
                                                                          'cmap': 'BrBG_r'},
    'timeMonthly_avg_activeTracerNonLocalTendency_salinityNonLocalTendency': {'label': 'Salt Non-local Tend (PSU/s)',
                                                                              'cmap': 'BrBG_r'}
}

MAX_MONTHS = 60  # 5 years for each simulation case

all_possible_sites = [
    "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
    "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
    "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4", "FSW1"
]


def get_field_data(ds_site, param, valid_layers_mask):
    if param == 'flow_speed':
        v_mer = 'timeMonthly_avg_velocityMeridional'
        v_zon = 'timeMonthly_avg_velocityZonal'
        if v_mer in ds_site and v_zon in ds_site:
            m_matrix = ds_site[v_mer].transpose('nVertLevels', 'Time').values[valid_layers_mask, :]
            z_matrix = ds_site[v_zon].transpose('nVertLevels', 'Time').values[valid_layers_mask, :]
            return np.sqrt(m_matrix ** 2 + z_matrix ** 2) * 100
    else:
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
                # Add offset based on actual months processed in Spin1
                res_patch_years.append((actual_spin1_months + spin6_months_cumulative) / 12.0)

    if not spin6_ds_list:
        print(f"⚠️ Warning: Missing Spin6 files for {res}. Plot will only show Spin1.")
        combined_datasets[res] = ds_spin1_cut
        patch_times_dict[res] = res_patch_years
        continue

    ds_spin6_full = xr.concat(spin6_ds_list, dim='Time')
    ds_spin6_cut = ds_spin6_full.isel(Time=slice(0, min(ds_spin6_full.sizes['Time'], MAX_MONTHS)))

    # --- Concatenate Along Time Dimension ---
    # We assign clean sequential data coordinates to the combined time dimension
    ds_combined = xr.concat([ds_spin1_cut, ds_spin6_cut], dim='Time', data_vars='minimal', coords='minimal',
                            compat='override')
    combined_datasets[res] = ds_combined

    # Store dynamic data patch structures
    patch_times_dict[res] = res_patch_years

# =============================================================================
# MAIN MULTI-DIMENSIONAL COORD LOOP
# =============================================================================
for site_name in all_possible_sites:

    site_exists_anywhere = False
    for res in resolutions:
        if res in combined_datasets and site_name in combined_datasets[res]['site'].values:
            site_exists_anywhere = True
            break

    if not site_exists_anywhere:
        continue

    for group_id, group_params in plotting_groups.items():
        num_cols = len(group_params)
        print(f" -> Site: {site_name} | Processing Group: {group_id}...")

        # Global scale normalization pass across concatenated timelines
        site_limits = {param: {'min': np.inf, 'max': -np.inf} for param in group_params}
        for res in resolutions:
            if res not in combined_datasets or site_name not in combined_datasets[res]['site'].values:
                continue

            ds_site_check = combined_datasets[res].sel(site=site_name)
            thick_t0 = ds_site_check['timeMonthly_avg_layerThickness'].isel(
                Time=0).values if 'timeMonthly_avg_layerThickness' in ds_site_check else None
            valid_mask = ~np.isnan(thick_t0) if thick_t0 is not None else slice(None)

            for param in group_params:
                vals = get_field_data(ds_site_check, param, valid_mask)
                if vals is not None and vals.size > 0 and not np.all(np.isnan(vals)):
                    site_limits[param]['min'] = min(site_limits[param]['min'], np.nanmin(vals))
                    site_limits[param]['max'] = max(site_limits[param]['max'], np.nanmax(vals))

        for param in group_params:
            if site_limits[param]['min'] == np.inf or np.isnan(site_limits[param]['min']):
                site_limits[param]['min'], site_limits[param]['max'] = 0.0, 1.0

        fig, axes = plt.subplots(nrows=len(resolutions), ncols=num_cols, figsize=(7 * num_cols, 16), sharex='col')
        if num_cols == 1:
            axes = axes[:, np.newaxis]

        fig.suptitle(f'Resolution Comparison Hovmöller Profiles — Site: {site_name} (Spin1 & Spin6 Combined)',
                     fontsize=18, fontweight='normal', y=0.96)

        for r_idx, res in enumerate(resolutions):
            axes[r_idx, 0].annotate(f'{res}', xy=(0, 0.5), xytext=(-55, 0),
                                    xycoords='axes fraction', textcoords='offset points',
                                    size=14, weight='normal', ha='right', va='center')

            if res not in combined_datasets:
                for c_idx in range(num_cols):
                    axes[r_idx, c_idx].text(0.5, 0.5, f"No Data available for {res}", ha='center', va='center')
                continue

            ds_res = combined_datasets[res]
            if site_name not in ds_res['site'].values:
                for c_idx in range(num_cols):
                    axes[r_idx, c_idx].text(0.5, 0.5, f"Site missing in {res}", ha='center', va='center')
                continue

            ds_site = ds_res.sel(site=site_name)
            total_months = ds_site.sizes['Time']

            # Build Continuous Combined Time Mapping Array
            time_coords, tick_positions, tick_labels = [], [], []
            for i in range(total_months):
                decimal_year = i / 12.0
                time_coords.append(decimal_year)

                # Dynamic labels generated cleanly every 12 months (Years 0-10)
                if i % 12 == 0:
                    tick_positions.append(decimal_year)
                    if i < MAX_MONTHS:
                        tick_labels.append(f"S1 Yr {i // 12 + 1}")
                    else:
                        tick_labels.append(f"S6 Yr {(i - MAX_MONTHS) // 12 + 1}")

            time_coords = np.array(time_coords)
            dx = 1.0 / 12.0
            x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

            # Build Vertical Grid
            ssh_t0 = float(
                ds_site['timeMonthly_avg_ssh'].isel(Time=0).values) if 'timeMonthly_avg_ssh' in ds_site else 0.0
            thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values
            valid_layers_mask = ~np.isnan(thick_t0)
            thick_valid = thick_t0[valid_layers_mask]

            if len(thick_valid) == 0:
                for c_idx in range(num_cols):
                    axes[r_idx, c_idx].text(0.5, 0.5, "No active ocean layers", ha='center', va='center')
                continue

            z_edges = ssh_t0 - np.insert(np.cumsum(thick_valid), 0, 0.0)

            # Plot each variable panel in the current row
            for c_idx, param in enumerate(group_params):
                ax = axes[r_idx, c_idx]
                meta = param_meta.get(param, {'label': param, 'cmap': 'viridis'})

                vmin = site_limits[param]['min']
                vmax = site_limits[param]['max']

                data_matrix = get_field_data(ds_site, param, valid_layers_mask)
                if data_matrix is None:
                    ax.text(0.5, 0.5, "Missing Fields", ha='center', va='center', color='red')
                    continue

                mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'], vmin=vmin, vmax=vmax,
                                     edgecolors='none')
                ax.grid(True, linestyle=':', alpha=0.3, color='black')

                # Solid red vertical line dividing Spin1 and Spin6 at precisely 5 years (60 months)
                if total_months > MAX_MONTHS:
                    ax.axvline(x=MAX_MONTHS / 12.0, color='magenta', linestyle='-.', linewidth=1.0, zorder=4,
                               label='Spin6 Start')

                # Minor internal netCDF segment break lines (if any exist inside the runs)
                for patch_time in patch_times_dict.get(res, []):
                    ax.axvline(x=patch_time, color='black', linestyle='--', linewidth=1.0, zorder=3)

                if r_idx == 0:
                    ax.set_title(meta['label'], fontsize=14, pad=12)
                if c_idx == 0:
                    ax.set_ylabel('Elevation z (m)', fontsize=11)

                if r_idx == len(resolutions) - 1:
                    ax.set_xticks(tick_positions)
                    ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=9)
                    ax.set_xlabel('Simulation Run Timeline', fontsize=12, labelpad=10)

                ax.set_xlim(x_edges[0], x_edges[-1])
                cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, aspect=18)
                cbar.ax.tick_params(labelsize=9)

        plt.subplots_adjust(hspace=0.22, wspace=0.18, bottom=0.12, top=0.89, left=0.12, right=0.95)
        out_img_name = os.path.join(output_dir, f'{site_name}_{group_id}_10Y.png')
        plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
        plt.close()

for res, open_ds in combined_datasets.items():
    open_ds.close()

print("\n Combined Spin1 & Spin6 Hovmöller profile processing completed successfully!")