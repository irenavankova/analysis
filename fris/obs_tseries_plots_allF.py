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

    # Temperature Tendencies (Using diverging colormaps centered around zero balance)
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

    # Salinity Tendencies (Using diverging colormaps centered around zero balance)
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

# Max length to plot: 5 years * 12 months = 60 months
MAX_MONTHS = 60

all_possible_sites = [
    "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
    "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
    "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4", "FSW1"
]


# Helper function to extract or compute the requested data array matrix
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
            # Check dimensions to safely transpose
            dim_order = ('nVertLevels', 'Time') if 'nVertLevels' in ds_site[param].dims else (ds_site[param].dims[0],
                                                                                              'Time')
            return ds_site[param].transpose(*dim_order).values[valid_layers_mask, :]
    return None


# =============================================================================
# MAIN MULTI-DIMENSIONAL COORD LOOP
# =============================================================================
for sim_case, res_dict in simulation_cases.items():
    print("\n" + "=" * 85)
    print(f"LOADING DATASETS FOR SIMULATION CASE: {sim_case}")
    print("=" * 85)

    case_datasets = {}
    patch_times_dict = {}

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

                cumulative_months += ds_part.sizes['Time']
                if idx < len(file_list) - 1:
                    patch_months.append(cumulative_months)

        if not datasets_to_combine:
            print(f"⚠️ Warning: Missing files for {res}_{sim_case}. Row will be skipped.")
            continue

        ds_full = xr.concat(datasets_to_combine, dim='Time')
        plot_months = min(ds_full.sizes['Time'], MAX_MONTHS)
        case_datasets[res] = ds_full.isel(Time=slice(0, plot_months))

        res_patch_years = []
        for m in patch_months:
            if m < plot_months:
                current_year = m // 12
                current_month_0indexed = m % 12
                decimal_year = current_year + (current_month_0indexed / 12.0)
                res_patch_years.append(decimal_year)

        patch_times_dict[res] = res_patch_years

    # -------------------------------------------------------------------------
    # ITERATE SITE BY SITE
    # -------------------------------------------------------------------------
    for site_name in all_possible_sites:

        site_exists_anywhere = False
        for res in resolutions:
            if res in case_datasets and site_name in case_datasets[res]['site'].values:
                site_exists_anywhere = True
                break

        if not site_exists_anywhere:
            continue

        # ---------------------------------------------------------------------
        # LOOP THROUGH EACH PLOTTING GROUP INDEPENDENTLY
        # ---------------------------------------------------------------------
        for group_id, group_params in plotting_groups.items():
            num_cols = len(group_params)
            print(f" -> Site: {site_name} | Processing Plot Grid: {group_id} ({num_cols} Columns)...")

            # Calculate uniform local bounds for each parameter in the current group
            site_limits = {param: {'min': np.inf, 'max': -np.inf} for param in group_params}

            for res in resolutions:
                if res not in case_datasets or site_name not in case_datasets[res]['site'].values:
                    continue

                ds_site_check = case_datasets[res].sel(site=site_name)
                thick_t0 = ds_site_check['timeMonthly_avg_layerThickness'].isel(
                    Time=0).values if 'timeMonthly_avg_layerThickness' in ds_site_check else None
                valid_mask = ~np.isnan(thick_t0) if thick_t0 is not None else slice(None)

                for param in group_params:
                    vals = get_field_data(ds_site_check, param, valid_mask)
                    if vals is not None and vals.size > 0 and not np.all(np.isnan(vals)):
                        site_limits[param]['min'] = min(site_limits[param]['min'], np.nanmin(vals))
                        site_limits[param]['max'] = max(site_limits[param]['max'], np.nanmax(vals))

            # Handle fallbacks safely
            for param in group_params:
                if site_limits[param]['min'] == np.inf or np.isnan(site_limits[param]['min']):
                    site_limits[param]['min'], site_limits[param]['max'] = 0.0, 1.0

            # Initialize custom dynamic layout canvas matrix: 4 rows x N parameters
            fig, axes = plt.subplots(nrows=len(resolutions), ncols=num_cols, figsize=(7 * num_cols, 16), sharex='col')

            # Re-wrap 1D axes arrays if a group only contains 1 variable to safeguard indexing logic
            if num_cols == 1:
                axes = axes[:, np.newaxis]

            fig.suptitle(
                f'Resolution Comparison Hovmöller Profiles — Site: {site_name} ({sim_case} - {group_id.upper()})',
                fontsize=18, fontweight='normal', y=0.96)

            # Build out canvas blocks
            for r_idx, res in enumerate(resolutions):

                # Handle row metrics string tagging on y-axis border
                axes[r_idx, 0].annotate(f'{res}', xy=(0, 0.5), xytext=(-55, 0),
                                        xycoords='axes fraction', textcoords='offset points',
                                        size=14, weight='normal', ha='right', va='center', rotation=0)

                if res not in case_datasets:
                    for c_idx in range(num_cols):
                        axes[r_idx, c_idx].text(0.5, 0.5, f"No Data available for {res}", ha='center', va='center')
                    continue

                ds_res = case_datasets[res]

                if site_name not in ds_res['site'].values:
                    for c_idx in range(num_cols):
                        axes[r_idx, c_idx].text(0.5, 0.5, f"Site missing in {res}", ha='center', va='center')
                    continue

                ds_site = ds_res.sel(site=site_name)
                plot_months = ds_site.sizes['Time']

                # Build Time / X-Axis Space
                time_coords, tick_positions, tick_labels = [], [], []
                for i in range(plot_months):
                    current_year = i // 12
                    current_month_0indexed = i % 12
                    decimal_year = current_year + (current_month_0indexed / 12.0)
                    time_coords.append(decimal_year)
                    if current_month_0indexed == 0:
                        tick_positions.append(decimal_year)
                        tick_labels.append(f"Yr {current_year} M1")

                time_coords = np.array(time_coords)
                dx = time_coords[1] - time_coords[0] if plot_months > 1 else 1.0 / 12.0
                x_edges = np.append(time_coords - dx / 2.0, time_coords[-1] + dx / 2.0)

                # Build Vertical Grid / Y-Axis Space
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

                    mesh = ax.pcolormesh(x_edges, z_edges, data_matrix, cmap=meta['cmap'],
                                         vmin=vmin, vmax=vmax, edgecolors='none')
                    ax.grid(True, linestyle=':', alpha=0.3, color='black')

                    for patch_time in patch_times_dict.get(res, []):
                        ax.axvline(x=patch_time, color='black', linestyle='-', linewidth=1.5, zorder=3)

                    if r_idx == 0:
                        ax.set_title(meta['label'], fontsize=14, pad=12, fontweight='normal')
                    if c_idx == 0:
                        ax.set_ylabel('Elevation z (m)', fontsize=11)

                    if r_idx == len(resolutions) - 1:
                        ax.set_xticks(tick_positions)
                        ax.set_xticklabels(tick_labels, rotation=0, ha='center', fontsize=10)
                        ax.set_xlabel('Model Year', fontsize=12, labelpad=10)

                    ax.set_xlim(x_edges[0], x_edges[-1])

                    cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, aspect=18)
                    cbar.ax.tick_params(labelsize=9)

            plt.subplots_adjust(hspace=0.20, wspace=0.18, bottom=0.08, top=0.89, left=0.12, right=0.95)
            out_img_name = os.path.join(output_dir, f'{site_name}_{group_id}_{sim_case}.png')
            plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
            plt.close()

    for res, open_ds in case_datasets.items():
        open_ds.close()

print("\n Combined resolution/parameter grid comparison plots completed successfully for all groups!")