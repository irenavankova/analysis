#!/usr/bin/env python3

import os
import glob
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# =====================================================================
# CONFIGURATION OPTIONS
# =====================================================================
DIRECTORY = "/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/LANL/obs4E3SM/subshelf_ctd/nc_files"
OUTPUT_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_plots/compare_obs"

# Model Folder and File Definitions
MODEL_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/obs"
MODEL_FILES = {
    "F1": "obs_tseries_F1_Spin6p1.nc",
    "F2": "obs_tseries_F2_Spin6p1.nc",
    "F4": "obs_tseries_F4_Spin6p1.nc",
    "F8": "obs_tseries_F8_Spin6p1.nc"
}

# Model Variable Target Configuration Strings
MODEL_T_VAR = "timeMonthly_avg_activeTracers_temperature"
MODEL_S_VAR = "timeMonthly_avg_activeTracers_salinity"
MODEL_THICK_VAR = "timeMonthly_avg_layerThickness"
MODEL_SSH_VAR = "timeMonthly_avg_ssh"

# Filter sites: Leave empty [] to plot ALL sites, or specify names to subset

INCLUDED_SITES = [
   "FNE1", "FNE3", "FSE1", "FSW1", "FSW2", "Site2", "Site3", "Site5"
]
# Save settings
SAVE_FIGURE = True
OUTPUT_DPI = 300

# ---------------------------------------------------------------------
# PLOTTING STYLES Configuration (Restored to Original Variant Colors)
# ---------------------------------------------------------------------
OBS_T_STYLE = {"color": "yellowgreen", "linestyle": "-", "linewidth": 2.0, "label": "Obs Temp"}
OBS_S_STYLE = {"color": "yellowgreen", "linestyle": "-", "linewidth": 2.0, "label": "Obs Sal"}

MODEL_STYLES = {
    "F1": {"color": "black", "linestyle": "-", "linewidth": 2.0, "label": "F1"},
    "F2": {"color": "orange", "linestyle": "-", "linewidth": 2.0, "label": "F2"},
    "F4": {"color": "dodgerblue", "linestyle": "-", "linewidth": 2.0, "label": "F4"},
    "F8": {"color": "brown", "linestyle": "-", "linewidth": 2.0, "label": "F8"}
}

MODEL_ALPHA = 0.08  # Shading transparency for individual model time points
MODEL_MEAN_LW = 2.0  # Line width for model time-averaged profiles

# Proportional aspect layout parameters (tall format)
FIG_WIDTH_PER_COL = 3.0
FIG_HEIGHT_PER_ROW = 4.0
MAX_COLUMNS = 3


# =====================================================================

def compute_model_z(ssh_series, thick_series):
    """
    Computes the z-level coordinates (mid-layer heights) for every time step.
    """
    n_time, n_levels = thick_series.shape
    z_model = np.zeros((n_time, n_levels))

    for t in range(n_time):
        ssh = ssh_series[t]
        thick = thick_series[t, :]
        interfaces = ssh - np.cumsum(thick)
        z_model[t, :] = interfaces + (thick / 2.0)

    return z_model


def main():
    # Find all observational profile netCDF files
    all_files = sorted(glob.glob(os.path.join(DIRECTORY, "*.nc")))
    nc_files = []

    for file_path in all_files:
        filename = os.path.basename(file_path)
        if "_" in filename:
            site_name = filename.split("_")[-1].replace(".nc", "")
        else:
            continue

        if not INCLUDED_SITES or site_name in INCLUDED_SITES:
            nc_files.append((file_path, site_name))

    if not nc_files:
        print(f"No matching observational netCDF files found in '{DIRECTORY}'.")
        return

    # Pre-load all model datasets into memory
    model_datasets = {}
    for key, filename in MODEL_FILES.items():
        path = os.path.join(MODEL_DIRECTORY, filename)
        if os.path.exists(path):
            model_datasets[key] = xr.open_dataset(path, decode_times=False)
        else:
            print(f"Warning: Model file {path} not found. Skipping {key}.")

    num_files = len(nc_files)
    cols = min(MAX_COLUMNS, num_files)
    rows = math.ceil(num_files / cols)

    # Helper sizing tool function for window framing engines
    def apply_screen_fit(fig_obj):
        try:
            backend = plt.get_backend()
            mgr = plt.get_current_fig_manager()
            if hasattr(mgr, 'window'):
                screen_width, screen_height = 14, 8
                if backend == 'TkAgg':
                    screen_width = mgr.window.winfo_screenwidth() / fig_obj.dpi
                    screen_height = mgr.window.winfo_screenheight() / fig_obj.dpi
                elif backend == 'MacOSX':
                    import AppKit
                    frame = AppKit.NSScreen.mainScreen().frame()
                    screen_width = frame.size.width / fig_obj.dpi
                    screen_height = frame.size.height / fig_obj.dpi

                aspect_ratio = (cols * FIG_WIDTH_PER_COL) / (rows * FIG_HEIGHT_PER_ROW)
                target_height = screen_height * 0.82
                target_width = target_height * aspect_ratio
                if target_width > screen_width:
                    target_width = screen_width * 0.82
                    target_height = target_width / aspect_ratio
                fig_obj.set_size_inches(target_width, target_height)
        except Exception:
            fig_obj.set_size_inches(9, 11)

    # =====================================================================
    # FIGURE 1: SALINITY PROFILES
    # =====================================================================
    fig_s, axes_s = plt.subplots(rows, cols, figsize=(cols * FIG_WIDTH_PER_COL, rows * FIG_HEIGHT_PER_ROW), sharey=True,
                                 constrained_layout=True)
    axes_s = np.array([axes_s]) if num_files == 1 else axes_s.flatten()
    legend_s = {}

    for i, (file_path, site_name) in enumerate(nc_files):
        ax = axes_s[i]

        # Plot Obs Salinity
        with xr.open_dataset(file_path, decode_times=False) as ds_obs:
            z_obs = np.squeeze(ds_obs['z'].values)
            s_obs = np.squeeze(ds_obs['salinity'].values)
            h_obs = ax.plot(s_obs, z_obs, **OBS_S_STYLE)[0]
            if i == 0: legend_s["Obs Sal"] = h_obs

        # Plot Model Salinity
        for key, ds_model in model_datasets.items():
            model_sites = [s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip() for s in
                           ds_model['site'].values]
            if site_name in model_sites:
                site_idx = model_sites.index(site_name)
                style = MODEL_STYLES[key]

                if MODEL_SSH_VAR in ds_model and MODEL_THICK_VAR in ds_model and MODEL_S_VAR in ds_model:
                    z_model = compute_model_z(ds_model[MODEL_SSH_VAR][:, site_idx].values,
                                              ds_model[MODEL_THICK_VAR][:, site_idx, :].values)
                    s_profiles = ds_model[MODEL_S_VAR][:, site_idx, :].values

                    for t in range(s_profiles.shape[0]):
                        ax.plot(s_profiles[t, :], z_model[t, :], color=style["color"], alpha=MODEL_ALPHA, linewidth=0.5)

                    h_m = ax.plot(np.nanmean(s_profiles, axis=0), np.nanmean(z_model, axis=0), color=style["color"],
                                  linestyle=style["linestyle"], linewidth=MODEL_MEAN_LW)[0]
                    if i == 0: legend_s[f"{key} Model"] = h_m

        ax.grid(True, linestyle=":", alpha=0.5)
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Salinity (PSU)")
        ax.tick_params(axis='x', labelbottom=True)  # Restored tick display explicitly

        # Handle label isolation rules
        if not (i >= (rows - 1) * cols or (i + cols) >= num_files):
            ax.set_xlabel("")

        ax.text(0.03, 0.1, f"{chr(97 + i)}) {site_name}", transform=ax.transAxes,
                fontsize=plt.rcParams['axes.labelsize'], weight='normal', va='top', ha='left', zorder=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.0))

    for j in range(i + 1, len(axes_s)): fig_s.delaxes(axes_s[j])
    fig_s.legend(legend_s.values(), legend_s.keys(), loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=5,
                 frameon=True)
    apply_screen_fit(fig_s)
    if SAVE_FIGURE: plt.savefig(os.path.join(OUTPUT_DIRECTORY, "comp_obs_profiles_Salinity.png"), dpi=OUTPUT_DPI,
                                bbox_inches='tight')

    # =====================================================================
    # FIGURE 2: TEMPERATURE PROFILES
    # =====================================================================
    fig_t, axes_t = plt.subplots(rows, cols, figsize=(cols * FIG_WIDTH_PER_COL, rows * FIG_HEIGHT_PER_ROW), sharey=True,
                                 constrained_layout=True)
    axes_t = np.array([axes_t]) if num_files == 1 else axes_t.flatten()
    legend_t = {}

    for i, (file_path, site_name) in enumerate(nc_files):
        ax = axes_t[i]

        # Plot Obs Temperature
        with xr.open_dataset(file_path, decode_times=False) as ds_obs:
            z_obs = np.squeeze(ds_obs['z'].values)
            t_obs = np.squeeze(ds_obs['potentialTemperature'].values)
            h_obs = ax.plot(t_obs, z_obs, **OBS_T_STYLE)[0]
            if i == 0: legend_t["Obs Temp"] = h_obs

        # Plot Model Temperature
        for key, ds_model in model_datasets.items():
            model_sites = [s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip() for s in
                           ds_model['site'].values]
            if site_name in model_sites:
                site_idx = model_sites.index(site_name)
                style = MODEL_STYLES[key]

                if MODEL_SSH_VAR in ds_model and MODEL_THICK_VAR in ds_model and MODEL_T_VAR in ds_model:
                    z_model = compute_model_z(ds_model[MODEL_SSH_VAR][:, site_idx].values,
                                              ds_model[MODEL_THICK_VAR][:, site_idx, :].values)
                    t_profiles = ds_model[MODEL_T_VAR][:, site_idx, :].values

                    for t in range(t_profiles.shape[0]):
                        ax.plot(t_profiles[t, :], z_model[t, :], color=style["color"], alpha=MODEL_ALPHA, linewidth=0.5)

                    h_m = ax.plot(np.nanmean(t_profiles, axis=0), np.nanmean(z_model, axis=0), color=style["color"],
                                  linestyle=style["linestyle"], linewidth=MODEL_MEAN_LW)[0]
                    if i == 0: legend_t[f"{key}"] = h_m

        ax.grid(True, linestyle=":", alpha=0.5)
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Potential Temperature (°C)")
        ax.tick_params(axis='x', labelbottom=True)  # Restored tick display explicitly

        if not (i >= (rows - 1) * cols or (i + cols) >= num_files):
            ax.set_xlabel("")

        ax.text(0.03, 0.1, f"{chr(97 + i)}) {site_name}", transform=ax.transAxes,
                fontsize=plt.rcParams['axes.labelsize'], weight='normal', va='top', ha='left', zorder=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.0))

    for j in range(i + 1, len(axes_t)): fig_t.delaxes(axes_t[j])
    fig_t.legend(legend_t.values(), legend_t.keys(), loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=5,
                 frameon=True)
    apply_screen_fit(fig_t)
    if SAVE_FIGURE: plt.savefig(os.path.join(OUTPUT_DIRECTORY, "comp_obs_profiles_Temperature.png"), dpi=OUTPUT_DPI,
                                bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()