#!/usr/bin/env python3

import os
import glob
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cftime
import matplotlib.ticker as ticker

# =====================================================================
# CONFIGURATION OPTIONS
# =====================================================================
DIRECTORY = "/Users/ivankova/Library/CloudStorage/GoogleDrive-irena.vanek@gmail.com/My Drive/Research/LANL/FRIS/development/transects_nc_files/FRIS/ApRES_timeseries"
OUTPUT_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_plots/compare_obs"

# Model Folder and File Definitions
MODEL_DIRECTORY = "/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/pts_tseries/obs"
MODEL_FILES = {
    "F1": "obs_tseries_F1_Spin6p1.nc",
    "F2": "obs_tseries_F2_Spin6p1.nc",
    "F4": "obs_tseries_F4_Spin6p1.nc",
    "F8": "obs_tseries_F8_Spin6p1.nc"
}

# CHOOSE MODEL VARIABLE TO PLOT:
# Option 1: "timeMonthly_avg_landIceFreshwaterFluxTotal"
# Option 2: "timeMonthly_avg_landIceFreshwaterFlux"
MODEL_VARIABLE = "timeMonthly_avg_landIceFreshwaterFlux"

# Filter sites: Leave empty [] to plot ALL sites, or specify names to subset
INCLUDED_SITES = [
    "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R15", "Site5", "FNE3", "FSE1", "FSW2"
]

# Save settings
SAVE_FIGURE = True
OUTPUT_DPI = 300

# Raw / Continuous observational melt line settings
RAW_COLOR = "gray"
RAW_LINESTYLE = "-"
RAW_ALPHA = 0.5
RAW_LINEWIDTH = 1
RAW_LABEL = "Obs Raw Melt"

# Monthly observational melt line settings
MONTHLY_COLOR = "yellowgreen"
MONTHLY_ALPHA = 0.2
MONTHLY_LINESTYLE = "-"
MONTHLY_LINEWIDTH = 2
MONTHLY_LABEL = "Obs Monthly Melt"

# ---------------------------------------------------------------------
# MODEL PLOTTING STYLES Configuration
# ---------------------------------------------------------------------
MODEL_STYLES = {
    "F1": {"color": "black", "linestyle": "-", "linewidth": 2.0, "label": "F1"},
    "F2": {"color": "orange", "linestyle": "-", "linewidth": 2.0, "label": "F2"},
    "F4": {"color": "dodgerblue", "linestyle": "-", "linewidth": 2.0, "label": "F4"},
    "F8": {"color": "brown", "linestyle": "-", "linewidth": 2.0, "label": "F8"}
}

# Figure settings
FIG_WIDTH_PER_COL = 4
FIG_HEIGHT_PER_ROW = 6
MAX_COLUMNS = 2  # Restructured to two columns only


# =====================================================================

def get_years_and_filter(time_var, melt_var):
    """
    Decodes the raw time variable allowing year zero, finds the dataset's end point,
    looks back 6 years, aligns the start boundary back to January 1st of that lookback year,
    filters the datasets, and outputs relative integer-aligned years.
    """
    days_raw = np.squeeze(time_var.values)
    melt_raw = np.squeeze(melt_var.values)

    units = time_var.attrs.get('units', 'days since 0000-01-01')
    dates = cftime.num2date(days_raw, units=units, calendar='standard', has_year_zero=True)

    end_date = dates[-1]
    target_start_year = end_date.year - 6

    jan_first_cftime = cftime.datetime(target_start_year, 1, 1, 0, 0, 0, calendar='standard', has_year_zero=True)
    jan_first_days = cftime.date2num(jan_first_cftime, units=units, calendar='standard', has_year_zero=True)

    mask = days_raw >= jan_first_days
    days_filtered = days_raw[mask]
    melt_filtered = melt_raw[mask]

    relative_years = (days_filtered - jan_first_days) / 365.25

    return relative_years, melt_filtered


def main():
    # Find all observational netCDF files
    all_files = sorted(glob.glob(os.path.join(DIRECTORY, "FRIS_basal_melt_*.nc")))
    nc_files = []

    for file_path in all_files:
        filename = os.path.basename(file_path)
        site_name = filename.replace("FRIS_basal_melt_", "").replace(".nc", "")

        if not INCLUDED_SITES or site_name in INCLUDED_SITES:
            nc_files.append((file_path, site_name))

    if not nc_files:
        print(f"No matching netCDF files found in '{DIRECTORY}' for the selected criteria.")
        return

    # Pre-load all model data datasets into memory
    model_datasets = {}
    for key, filename in MODEL_FILES.items():
        path = os.path.join(MODEL_DIRECTORY, filename)
        if os.path.exists(path):
            model_datasets[key] = xr.open_dataset(path)
        else:
            print(f"Warning: Model file {path} not found. Skipping {key}.")

    num_files = len(nc_files)
    cols = min(MAX_COLUMNS, num_files)
    rows = math.ceil(num_files / cols)

    fig, axes = plt.subplots(
        rows, cols,
        figsize=(cols * FIG_WIDTH_PER_COL, rows * FIG_HEIGHT_PER_ROW),
        sharex=True, sharey=False,
        constrained_layout=True
    )

    if num_files == 1:
        axes = np.array([axes])
    else:
        axes = axes.flatten()

    for i, (file_path, site_name) in enumerate(nc_files):
        ax = axes[i]

        # 1. Plot Observational Timeseries Data
        with xr.open_dataset(file_path, decode_times=False) as ds:
            years_raw, melt_raw = get_years_and_filter(ds['time'], ds['melt_timeseries'])
            years_monthly, melt_monthly = get_years_and_filter(ds['time_monthly'], ds['melt_timeseries_monthly'])

            ax.plot(
                years_monthly, melt_monthly,
                color=MONTHLY_COLOR, linestyle=MONTHLY_LINESTYLE, linewidth=MONTHLY_LINEWIDTH,
                label=MONTHLY_LABEL if i == 0 else ""
            )

        # 2. Plot Model Outputs if available for this specific site name
        for key, ds_model in model_datasets.items():
            model_sites = [s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip()
                           for s in ds_model['site'].values]

            if site_name in model_sites:
                site_idx = model_sites.index(site_name)

                if MODEL_VARIABLE in ds_model:
                    flux_data = ds_model[MODEL_VARIABLE][:, site_idx].values
                else:
                    print(f"Warning: Variable '{MODEL_VARIABLE}' not found in {key} dataset.")
                    continue

                melt_model_ma = (flux_data / 917.0) * 31536000.0
                model_time_years = np.linspace(0, 6, len(flux_data))

                style = MODEL_STYLES[key]
                ax.plot(
                    model_time_years, melt_model_ma,
                    color=style["color"], linestyle=style["linestyle"], linewidth=style["linewidth"],
                    label=style["label"] if i == 0 else ""
                )

        # Design cosmetics
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.set_ylabel("Melt Rate (m/a)")
        ax.set_xlim(0, 7)  # Extended x-axis range to 7 years

        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        # Only set x-label on bottom-most active subplots
        if i >= (rows - 1) * cols or (i + cols) >= num_files:
            ax.set_xlabel("Time (years from January 1st)")

        # Alphabetical index generation (a, b, c...)
        letter_prefix = chr(97 + i)  # ASCII 97 starts at 'a'
        combined_label = f"{letter_prefix}) {site_name}"

        # Inner label styled identically to default y-label parameters (normal weight)
        ax.text(
            0.86, 0.1, combined_label,
            transform=ax.transAxes,
            fontsize=plt.rcParams['axes.labelsize'],
            weight='normal',
            va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.0)
        )

    # Hide extra empty subplots in grid frame
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=6, frameon=True)

    # --- Screen Fit Engine Override ---
    # --- Screen Fit Engine Override ---
    try:
        backend = plt.get_backend()
        mgr = plt.get_current_fig_manager()

        if hasattr(mgr, 'window'):
            # 1. Determine your physical monitor's usable height
            screen_width, screen_height = 14, 8  # Default fallbacks

            if backend == 'TkAgg':
                screen_width = mgr.window.winfo_screenwidth() / fig.dpi
                screen_height = mgr.window.winfo_screenheight() / fig.dpi
            elif backend == 'MacOSX':
                import AppKit
                frame = AppKit.NSScreen.mainScreen().frame()
                screen_width = frame.size.width / fig.dpi
                screen_height = frame.size.height / fig.dpi
            elif 'Qt' in backend:
                try:
                    from PyQt5.QtWidgets import QDesktopWidget
                    geom = QDesktopWidget().screenGeometry()
                    screen_width = geom.width() / fig.dpi
                    screen_height = geom.height() / fig.dpi
                except ImportError:
                    pass

            # 2. Scale the figure so it fits the screen height but preserves your aspect ratio
            desired_width = cols * FIG_WIDTH_PER_COL
            desired_height = rows * FIG_HEIGHT_PER_ROW
            aspect_ratio = desired_width / desired_height

            # Use 85% of screen height so it doesn't get cut off by Mac menu/dock bars
            target_height = screen_height * 0.85
            target_width = target_height * aspect_ratio

            # If the calculated width is somehow wider than the screen, scale by width instead
            if target_width > screen_width:
                target_width = screen_width * 0.85
                target_height = target_width / aspect_ratio

            fig.set_size_inches(target_width, target_height)

            # If using native window managers, don't use 'showMaximized' as it forces stretching
            if backend == 'TkAgg':
                mgr.window.geometry(f"{int(target_width * fig.dpi)}x{int(target_height * fig.dpi)}+50+50")
        else:
            fig.set_size_inches(10, 12)
    except Exception:
        fig.set_size_inches(10, 12)


    if SAVE_FIGURE:
        short_var_name = MODEL_VARIABLE.replace("timeMonthly_avg_", "")
        dynamic_filename = f"comp_obs_apares_{short_var_name}.png"
        output_path = os.path.join(OUTPUT_DIRECTORY, dynamic_filename)

        plt.savefig(output_path, dpi=OUTPUT_DPI, bbox_inches='tight')
        print(f"Figure successfully saved to: {output_path}")

    plt.show()


if __name__ == "__main__":
    main()