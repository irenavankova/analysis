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

# Filter sites: Leave empty [] to plot ALL sites, or specify names to subset
INCLUDED_SITES = [
    "FNE1", "FNE3", "FSE1", "FSW1", "FSW2", "Site2", "Site3", "Site5"
]

# Save settings
SAVE_FIGURE = True
OUTPUT_DPI = 300

# ---------------------------------------------------------------------
# PLOTTING STYLES Configuration
# ---------------------------------------------------------------------
OBS_STYLE = {"color": "yellowgreen", "linestyle": "-", "linewidth": 2.5, "label": "Obs CTD"}

MODEL_STYLES = {
    "F1": {"color": "black", "linestyle": "-", "label": "F1"},
    "F2": {"color": "orange", "linestyle": "-", "label": "F2"},
    "F4": {"color": "dodgerblue", "linestyle": "-", "label": "F4"},
    "F8": {"color": "brown", "linestyle": "-", "label": "F8"}
}

MODEL_ALPHA = 0.08  # Shading transparency for individual model time points
MODEL_MEAN_LW = 2.0  # Line width for model time-averaged profiles

# Proportional aspect layout parameters (T-S diagrams look best with uniform square aspects)
FIG_WIDTH_PER_COL = 4.5
FIG_HEIGHT_PER_ROW = 4.5
MAX_COLUMNS = 3


# =====================================================================

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

    # Setup the single joint canvas
    fig, axes = plt.subplots(
        rows, cols,
        figsize=(cols * FIG_WIDTH_PER_COL, rows * FIG_HEIGHT_PER_ROW),
        sharex=False, sharey=False,
        constrained_layout=True
    )

    if num_files == 1:
        axes = np.array([axes])
    else:
        axes = axes.flatten()

    legend_handles = {}

    for i, (file_path, site_name) in enumerate(nc_files):
        ax = axes[i]

        # 1. Plot Observational Profiles in T-S Space
        with xr.open_dataset(file_path, decode_times=False) as ds_obs:
            t_obs = np.squeeze(ds_obs['potentialTemperature'].values)
            s_obs = np.squeeze(ds_obs['salinity'].values)

            h_obs = ax.plot(s_obs, t_obs, **OBS_STYLE)[0]
            if i == 0:
                legend_handles["Obs CTD"] = h_obs

        # 2. Plot Model Datasets in T-S Space
        for key, ds_model in model_datasets.items():
            model_sites = [s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip()
                           for s in ds_model['site'].values]

            if site_name in model_sites:
                site_idx = model_sites.index(site_name)
                style = MODEL_STYLES[key]

                if MODEL_T_VAR in ds_model and MODEL_S_VAR in ds_model:
                    t_profiles = ds_model[MODEL_T_VAR][:, site_idx, :].values  # Shape: (Time, nVertLevels)
                    s_profiles = ds_model[MODEL_S_VAR][:, site_idx, :].values  # Shape: (Time, nVertLevels)

                    # Plot each time frame profile line path across vertical layers
                    for t in range(t_profiles.shape[0]):
                        ax.plot(s_profiles[t, :], t_profiles[t, :], color=style["color"], alpha=MODEL_ALPHA,
                                linewidth=0.5)

                    # Calculate and plot the joint Temporal Mean Profile line
                    t_mean = np.nanmean(t_profiles, axis=0)
                    s_mean = np.nanmean(s_profiles, axis=0)

                    h_m = ax.plot(s_mean, t_mean, color=style["color"], linestyle=style["linestyle"],
                                  linewidth=MODEL_MEAN_LW)[0]
                    if i == 0:
                        legend_handles[f"Model {key}"] = h_m

        # Layout styling mechanics
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.set_xlabel("Salinity (PSU)")
        ax.set_ylabel("Potential Temperature (°C)")
        ax.tick_params(axis='both', labelbottom=True, labelleft=True)  # Global ticks configured manually

        # Suppress labels if outside border cells
        is_bottom_row = i >= (rows - 1) * cols or (i + cols) >= num_files
        is_first_col = (i % cols == 0)

        if not is_bottom_row:
            ax.set_xlabel("")
        if not is_first_col:
            ax.set_ylabel("")

        # Subplot text tracking labels
        ax.text(
            0.03, 0.95, f"{chr(97 + i)}) {site_name}",
            transform=ax.transAxes,
            fontsize=plt.rcParams['axes.labelsize'],
            weight='normal',
            va='top', ha='left', zorder=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.7)
        )

    # Remove trailing empty frame allocations
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Centered Legend header layout orientation Pinned above the ceiling
    fig.legend(legend_handles.values(), legend_handles.keys(), loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=5,
               frameon=True)

    # --- Screen Fit Engine Proportional Aspect Ratio Override ---
    try:
        backend = plt.get_backend()
        mgr = plt.get_current_fig_manager()
        if hasattr(mgr, 'window'):
            screen_width, screen_height = 14, 8
            if backend == 'TkAgg':
                screen_width = mgr.window.winfo_screenwidth() / fig.dpi
                screen_height = mgr.window.winfo_screenheight() / fig.dpi
            elif backend == 'MacOSX':
                import AppKit
                frame = AppKit.NSScreen.mainScreen().frame()
                screen_width = frame.size.width / fig.dpi
                screen_height = frame.size.height / fig.dpi

            aspect_ratio = (cols * FIG_WIDTH_PER_COL) / (rows * FIG_HEIGHT_PER_ROW)
            target_height = screen_height * 0.82
            target_width = target_height * aspect_ratio
            if target_width > screen_width:
                target_width = screen_width * 0.82
                target_height = target_width / aspect_ratio
            fig.set_size_inches(target_width, target_height)
    except Exception:
        fig.set_size_inches(11, 8)

    if SAVE_FIGURE:
        plt.savefig(os.path.join(OUTPUT_DIRECTORY, "comp_obs_TS_diagrams.png"), dpi=OUTPUT_DPI, bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()