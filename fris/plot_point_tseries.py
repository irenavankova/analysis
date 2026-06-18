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

ts_name = 'pts_berknerwest'

# >>> USER SWITCH: Set to 1 to save to disk, or 0 to display directly on screen <<<
OPT_SAVE = 1

input_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/{ts_name}'
output_dir = f'/Users/ivankova/Desktop/Fris_hr/Fris_plots/point_tseries/{ts_name}'
if OPT_SAVE == 1:
    os.makedirs(output_dir, exist_ok=True)

# Centralized constants for unit conversion
sec_per_year = 365.25 * 24 * 3600  # ~31,557,600 seconds
rho_fw = 1000.0  # Freshwater density (kg/m^3)
m2cm = 100

# Centralized line property definitions
Lwides1 = 2.0
Lwides6 = 1.5
lstyl1 = '--'
lstyl6 = '-'

S6_start_yr = 0

# Target dates for vertical lines (defined as tuple: (year, month))
vline_targets = [(2, 11), (0, 12), (4, 6)]
vline_colors = ['yellowgreen', 'darkgray', 'darkgray']

# --- Manual Axis Limits ---
vmin_temp, vmax_temp = -2.5, -1.3  # Min/Max for Temperature (°C)
vmin_salt, vmax_salt = 34.4, 35.0  # Min/Max for Salinity (PSU)

# >>> NEW CONFIGURATION: Separate mapping dictionary for styling properties <<<
style_config = {
    'Spin6': {
        'F8': {'color': 'brown',      'linewidth': Lwides6, 'linestyle': lstyl6},
        'F4': {'color': 'dodgerblue', 'linewidth': Lwides6, 'linestyle': lstyl6},
        'F2': {'color': 'orange',     'linewidth': Lwides6, 'linestyle': lstyl6},
        'F1': {'color': 'black',      'linewidth': Lwides6, 'linestyle': lstyl6},
    },
    'Spin1': {
        'F8': {'color': 'burlywood',          'linewidth': Lwides1, 'linestyle': lstyl1},
        'F4': {'color': 'lightskyblue', 'linewidth': Lwides1, 'linestyle': lstyl1},
        'F2': {'color': 'gold',       'linewidth': Lwides1, 'linestyle': lstyl1},
        'F1': {'color': 'lightgray',         'linewidth': Lwides1, 'linestyle': lstyl1},
    }
}

# Raw file list properties mapping tracking timeline baselines
files_raw = [
    {'filename': f'{ts_name}_tseries_F8_Spin1p1.nc', 'res': 'F8', 'spin': 'Spin1', 'label': 'F8 (Spin1)', 'start_year': 0, 'start_month': 1,  'in_legend': True},
    {'filename': f'{ts_name}_tseries_F8_Spin6p1.nc', 'res': 'F8', 'spin': 'Spin6', 'label': 'F8 (Spin6)',   'start_year': S6_start_yr, 'start_month': 1, 'in_legend': True},
    {'filename': f'{ts_name}_tseries_F4_Spin1p1.nc', 'res': 'F4', 'spin': 'Spin1', 'label': 'F4 (Spin1)', 'start_year': 0, 'start_month': 1,  'in_legend': True},
    {'filename': f'{ts_name}_tseries_F4_Spin6p1.nc', 'res': 'F4', 'spin': 'Spin6', 'label': 'F4 (Spin6)',   'start_year': S6_start_yr, 'start_month': 1, 'in_legend': True},
    {'filename': f'{ts_name}_tseries_F2_Spin1p1.nc', 'res': 'F2', 'spin': 'Spin1', 'label': 'F2 (Spin1)', 'start_year': 0, 'start_month': 1,  'in_legend': True},
    {'filename': f'{ts_name}_tseries_F2_Spin1p2.nc', 'res': 'F2', 'spin': 'Spin1', 'label': 'F2 (Spin1p2)', 'start_year': 2, 'start_month': 11, 'in_legend': False},
    {'filename': f'{ts_name}_tseries_F2_Spin6p1.nc', 'res': 'F2', 'spin': 'Spin6', 'label': 'F2 (Spin6)',   'start_year': S6_start_yr, 'start_month': 1, 'in_legend': True},
    {'filename': f'{ts_name}_tseries_F1_Spin1p1.nc', 'res': 'F1', 'spin': 'Spin1', 'label': 'F1 (Spin1)', 'start_year': 0, 'start_month': 1,  'in_legend': True},
    {'filename': f'{ts_name}_tseries_F1_Spin1p2.nc', 'res': 'F1', 'spin': 'Spin1', 'label': 'F1 (Spin1p2)', 'start_year': 0, 'start_month': 12, 'in_legend': False},
    {'filename': f'{ts_name}_tseries_F1_Spin1p3.nc', 'res': 'F1', 'spin': 'Spin1', 'label': 'F1 (Spin1p3)', 'start_year': 4, 'start_month': 6,  'in_legend': False},
    {'filename': f'{ts_name}_tseries_F1_Spin6p1.nc', 'res': 'F1', 'spin': 'Spin6', 'label': 'F1 (Spin6)',   'start_year': S6_start_yr, 'start_month': 1, 'in_legend': True},
]

# Build the complete files_config dynamically using style mappings
files_config = []
for f in files_raw:
    style = style_config[f['spin']][f['res']]
    f.update(style)
    files_config.append(f)

# Parameters Metadata
param_meta = {
    'timeMonthly_avg_landIceFreshwaterFlux': {
        'label': 'Melt rate (m/a)',
        'vmin': None,
        'vmax': None,
        'target_layer': 'none'
    },
    'timeMonthly_avg_landIceFrictionVelocity': {
        'label': 'Ustar (cm/s)',
        'vmin': None,
        'vmax': None,
        'target_layer': 'none'
    },
    'timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature': {
        'label': 'BL temperature (°C)',
        'vmin': None,
        'vmax': None,
        'target_layer': 'none'
    },
    'timeMonthly_avg_activeTracers_temperature': {
        'label': 'Temperature (°C)',
        'vmin': None,
        'vmax': None,
        'target_layer': 'bottom'
    },
    'timeMonthly_avg_activeTracers_salinity': {
        'label': 'Salinity (PSU)',
        'vmin': None,
        'vmax': None,
        'target_layer': 'bottom'
    }
}

all_possible_sites = []
if ts_name == 'obs':
    all_possible_sites = [
        "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
        "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
        "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4", "FSW1"
    ]

vline_positions = [yr + ((mo - 1) / 12.0) for yr, mo in vline_targets]


def get_targeted_layer_data(ds_site, param, valid_layers_mask, layer_target):
    """Extracts top, bottom, or un-dimensioned 2D variable data cleanly."""
    if param in ds_site:
        variable = ds_site[param]

        if 'nVertLevels' not in variable.dims and layer_target.lower() == 'none':
            return variable.values.copy()

        dim_order = ('nVertLevels', 'Time') if 'nVertLevels' in variable.dims else (variable.dims[0], 'Time')
        full_matrix = variable.transpose(*dim_order).values
        valid_matrix = full_matrix[valid_layers_mask, :]

        if valid_matrix.shape[0] > 0:
            if layer_target.lower() == 'top':
                return valid_matrix[0, :].copy()
            elif layer_target.lower() == 'bottom':
                return valid_matrix[-1, :].copy()
    return None


# =============================================================================
# DATA VALIDATION BLOCK
# =============================================================================
print("\nValidating explicitly configured files...")
valid_files = []
for cfg in files_config:
    full_path = os.path.join(input_dir, cfg['filename'])
    if os.path.exists(full_path):
        cfg['path'] = full_path
        valid_files.append(cfg)

        if not all_possible_sites:
            with xr.open_dataset(full_path) as ds_test:
                if 'site' in ds_test:
                    all_possible_sites = list(ds_test['site'].values)
                    print(f" Found {len(all_possible_sites)} sites dynamically.")
    else:
        print(f"Warning: File not found, skipping: {cfg['filename']}")

if not valid_files:
    raise FileNotFoundError("Error: None of the configured NetCDF files were found.")

# =============================================================================
# MAIN PLOTTING LOOP (Per Site)
# =============================================================================
for site_name in all_possible_sites:

    site_exists = False
    for cfg in valid_files:
        with xr.open_dataset(cfg['path']) as ds:
            if 'site' in ds and site_name in ds['site'].values:
                site_exists = True
                break
    if not site_exists:
        continue

    print(f" -> Generating Time-series Figure for Site: {site_name}...")

    fig, axes = plt.subplots(nrows=len(param_meta), ncols=1, figsize=(11, 3.5 * len(param_meta)), sharex=True)
    if len(param_meta) == 1:
        axes = [axes]

    fig.suptitle(f'Site: {site_name}',
                 fontsize=12, fontweight='normal')

    global_max_x = 0
    legend_handles = []
    legend_labels = []

    # Loop over parameters (subplots)
    for p_idx, (param, meta) in enumerate(param_meta.items()):
        ax = axes[p_idx]

        target_layer = meta.get('target_layer', 'top').lower()

        if target_layer == 'none':
            ax.set_ylabel(f"{meta['label']}", fontsize=10)
        else:
            ax.set_ylabel(f"{meta['label']}\n[{target_layer.upper()} Layer]", fontsize=10)

        if meta['vmin'] is not None and meta['vmax'] is not None:
            ax.set_ylim(meta['vmin'], meta['vmax'])

        ax.grid(True, linestyle=':', alpha=0.5, color='gray')

        #for x_pos, line_color in zip(vline_positions, vline_colors):
        #    ax.axvline(x=x_pos, color=line_color, linestyle=':', linewidth=1.2, alpha=0.8, zorder=1)

        for cfg in valid_files:
            with xr.open_dataset(cfg['path']) as ds:
                if 'site' not in ds or site_name not in ds['site'].values:
                    continue
                if param not in ds.data_vars:
                    continue

                ds_site = ds.sel(site=site_name)
                num_timesteps = ds_site.sizes['Time']

                start_year = cfg['start_year']
                start_month = cfg['start_month']
                time_coords = []
                for i in range(num_timesteps):
                    total_months = (start_month - 1) + i
                    current_year = start_year + (total_months // 12)
                    current_month_0indexed = total_months % 12
                    decimal_year = current_year + (current_month_0indexed / 12.0)
                    time_coords.append(decimal_year)

                time_coords = np.array(time_coords)
                global_max_x = max(global_max_x, time_coords[-1])

                if 'timeMonthly_avg_layerThickness' in ds_site:
                    thick_t0 = ds_site['timeMonthly_avg_layerThickness'].isel(Time=0).values
                    valid_layers_mask = ~np.isnan(thick_t0)
                else:
                    valid_layers_mask = None

                # Extract the 1D timeline array
                layer_data = get_targeted_layer_data(ds_site, param, valid_layers_mask, layer_target=target_layer)
                if layer_data is None:
                    continue

                # Apply conversions based on the current field parameter
                if param == 'timeMonthly_avg_landIceFreshwaterFlux':
                    layer_data = layer_data * sec_per_year / rho_fw

                if param == 'timeMonthly_avg_landIceFrictionVelocity':
                    layer_data = layer_data * m2cm  # Fixed from division to multiplication for m/s -> cm/s

                line, = ax.plot(time_coords, layer_data,
                                label=cfg['label'],
                                color=cfg['color'],
                                linewidth=cfg['linewidth'],
                                linestyle=cfg['linestyle'],
                                alpha=0.85)

                if cfg.get('in_legend', True) and cfg['label'] not in legend_labels:
                    legend_handles.append(line)
                    legend_labels.append(cfg['label'])

    bottom_ax = axes[-1]
    bottom_ax.set_xlabel("Adjusted Model Year", fontsize=11, labelpad=8)
    bottom_ax.set_xlim(0, global_max_x)

    if legend_handles:
        fig.legend(legend_handles, legend_labels, loc='upper center',
                   bbox_to_anchor=(0.5, 0.94), ncol=4, frameon=True, fontsize=10)

    plt.subplots_adjust(hspace=0.22, bottom=0.10, top=0.88, left=0.12, right=0.95)

    # >>> CRITICAL ACTION BLOCK: Handle saving vs direct image showing <<<
    if OPT_SAVE == 1:
        out_img_name = os.path.join(output_dir, f'{site_name}_full_mixed_dimensions_comparison.png')
        plt.savefig(out_img_name, dpi=200, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

print("\nProcessing workflow execution concluded successfully!")