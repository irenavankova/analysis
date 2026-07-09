#!/usr/bin/env python3
import os
import glob
import numpy as np
import xarray as xr
from multiprocessing import Pool

import matplotlib

matplotlib.use('Agg')  # Force non-interactive backend for cluster environments
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm  # Added for logarithmic histogram scaling
import gsw
import gmask_reg
import cmocean


# =========================================================================
# 1. Parallel Worker Function for a Single Simulation Configuration
# =========================================================================
def process_single_ts_task(args):
    """Processes TS diagrams for all specified regions and seasons for a single simulation resolution."""
    Fnum, cases, RUN_TYPE, TARGET_YEARS, regions_to_plot, dir_fig_save, TS_bg_config, seasonal_windows = args
    dx = f'F{Fnum}'

    # Resolve case string identifiers for filenames (Matches plot_spatial_stats.py logic)
    cases_processed = []
    for sec, subsec in cases:
        subsec_str = subsec if (sec == 'Spin1' or subsec != 'p1') else ''
        cases_processed.append(f"{sec}{subsec_str}")
    combined_cases_str = "_".join(cases_processed)

    print(f"=========================================================================\n"
          f" Starting Execution: {dx}_{combined_cases_str} | Regions: {regions_to_plot}\n"
          f"=========================================================================")

    # Path construction following plot_spatial_stats.py
    run_name_mask = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
    fpath_mask = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name_mask}/run'
    mesh_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

    if not os.path.exists(mesh_file):
        print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}. Skipping task.")
        return

    # Retrieve Region Masks from gmask_reg.py for this specific mesh
    iam = gmask_reg.get_mask(regions_to_plot, mesh_file, opt_noGL=0, opt_wct=1)

    # Reconstruct base directories to gather monthly netCDF outputs (From plot_spatial_stats.py)
    unique_months_dict = {}
    for sec, subsec in cases:
        if sec == 'Spin6':
            if Fnum == '8' and subsec == 'GMF1':
                run_name = "20240503.GMPAS-JRA1p5-DIB-PISMF-DGMHT.TL319_FRISwISC08to60E3r1.spinY6_GMF1.chicoma-cpu"
                fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'
            else:
                run_name = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
                fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'
        elif sec == 'Spin1':
            if Fnum == '8':
                run_name = "20231114.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC08to60E3r1.spinup.chicoma-cpu"
            elif Fnum == '4':
                run_name = "20231108.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC04to60E3r1.spinup.chicoma-cpu"
            elif Fnum == '2':
                run_name = "20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.chicoma-cpu" if subsec == 'p1' else "20231208.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.anvil"
            elif Fnum == '1':
                if subsec == 'p1':
                    run_name = "20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.chicoma-cpu"
                elif subsec == 'p2':
                    run_name = "20231209.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.anvil"
                else:
                    run_name = "20240201.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinupY5.chicoma-cpu"
            fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY1/{run_name}/run'

        for yr_str in TARGET_YEARS:
            try:
                yr_int = int(yr_str)
                file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.{yr_int:04d}-*-*.nc"
                found_files = sorted(glob.glob(file_pattern))
                for file_path in found_files:
                    date_part = os.path.basename(file_path).split('.')[-2]
                    unique_months_dict[date_part] = file_path
            except (IndexError, ValueError):
                continue

    year_file_list = [unique_months_dict[k] for k in sorted(unique_months_dict.keys())]

    if not year_file_list:
        print(f"--> Warning: No monthly files found matching target years {TARGET_YEARS} for {dx}. Skipping.")
        return

    # Load cell grid geometry parameters once per resolution
    with xr.open_dataset(mesh_file) as dsM:
        areaCell = dsM['areaCell'].values  # Shape: (nCells,)
        maxLevelCell = dsM['maxLevelCell'].values - 1  # Shape: (nCells,)

    # Unpack TS background parameters for mapping
    y_lim = TS_bg_config['y_lim']
    x_lim = TS_bg_config['x_lim']
    PSbins = TS_bg_config['PSbins']
    PSgrid = TS_bg_config['PSgrid']
    PTgrid = TS_bg_config['PTgrid']
    neutralDensity = TS_bg_config['neutralDensity']
    contours = TS_bg_config['contours']
    PTFreezing = TS_bg_config['PTFreezing']

    years_str = f"Years_{TARGET_YEARS[0]}-{TARGET_YEARS[-1]}" if len(TARGET_YEARS) > 1 else f"Year_{TARGET_YEARS[0]}"

    # =========================================================================
    # Loop over the requested temporal windows (Annual + Seasonal subsets)
    # =========================================================================
    for label, allowed_months in seasonal_windows.items():

        # Filter files belonging to targeted months
        filtered_file_list = []
        for file_path in year_file_list:
            date_part = os.path.basename(file_path).split('.')[-2]  # e.g., '0002-03-01_00000'
            try:
                month_val = int(date_part.split('-')[1])  # extracts standard MM index
                if allowed_months is None or month_val in allowed_months:
                    filtered_file_list.append(file_path)
            except (IndexError, ValueError):
                continue

        if not filtered_file_list:
            print(f"--> Warning: No matching season files for [{label}] under resolution {dx}. Skipping window.")
            continue

        # Load and temporally average 3D properties across selected subset of months
        temp_list, salt_list, thick_list = [], [], []
        for file_path in filtered_file_list:
            with xr.open_dataset(file_path) as ds:
                temp_list.append(ds['timeMonthly_avg_activeTracers_temperature'].isel(Time=0).values)
                salt_list.append(ds['timeMonthly_avg_activeTracers_salinity'].isel(Time=0).values)
                thick_list.append(ds['timeMonthly_avg_layerThickness'].isel(Time=0).values)

        PT_mean = np.mean(np.array(temp_list), axis=0)  # Shape: (nCells, nVertLevels)
        PS_mean = np.mean(np.array(salt_list), axis=0)  # Shape: (nCells, nVertLevels)
        H_mean = np.mean(np.array(thick_list), axis=0)  # Shape: (nCells, nVertLevels)

        # Process and generate a plot for each region sequentially inside this window
        for r_idx, region_name in enumerate(regions_to_plot):
            region_mask = iam[r_idx, :]  # Shape: (nCells,)

            if not np.any(region_mask):
                continue

            # Extract only horizontal columns belonging to the current region mask
            PT_reg = PT_mean[region_mask, :]  # Shape: (nCells_in_reg, nVertLevels)
            PS_reg = PS_mean[region_mask, :]  # Shape: (nCells_in_reg, nVertLevels)
            H_reg = H_mean[region_mask, :]  # Shape: (nCells_in_reg, nVertLevels)
            areaCell_reg = areaCell[region_mask]  # Shape: (nCells_in_reg,)
            maxLevel_reg = maxLevelCell[region_mask]  # Shape: (nCells_in_reg,)

            # Calculate exact 3D grid volumes using broadcasting
            volume_reg = H_reg * areaCell_reg[:, np.newaxis]

            # Construct custom 2D vertical mask bounded by maxLevelCell per region column
            num_cells_reg, num_levels = PT_reg.shape
            level_indices = np.arange(num_levels)[np.newaxis, :]  # Shape: (1, nVertLevels)
            valid_vertical_mask = level_indices <= maxLevel_reg[:, np.newaxis]  # Shape: (nCells_in_reg, nVertLevels)

            # Flatten arrays safely using the 2D vertical indices mask
            PT_flat = PT_reg[valid_vertical_mask]
            PS_flat = PS_reg[valid_vertical_mask]
            Vol_flat = volume_reg[valid_vertical_mask]

            # Clean NaN data values
            nan_mask = np.isnan(PT_flat) | np.isnan(PS_flat) | np.isnan(Vol_flat)
            PT_flat = PT_flat[~nan_mask]
            PS_flat = PS_flat[~nan_mask]
            Vol_flat = Vol_flat[~nan_mask]

            if len(PT_flat) == 0:
                continue

            # Volume-weighted core spatial mean calculation
            PT_core_avg = np.dot(PT_flat, Vol_flat) / np.sum(Vol_flat)
            PS_core_avg = np.dot(PS_flat, Vol_flat) / np.sum(Vol_flat)

            # -----------------------------------------------------------------
            # Render Figure (Explicit figure flushing to avoid memory leaks)
            # -----------------------------------------------------------------
            fig, ax = plt.subplots(figsize=(6, 5))  # Slightly widened to comfortably fit the colorbar

            # Plot background potential density contours
            CS = ax.contour(PSgrid, PTgrid, neutralDensity, contours, linestyles=':', linewidths=0.5, colors='k',
                            zorder=2)
            ax.clabel(CS, fontsize=8, inline=1, fmt='%4.2f')

            # Surface Freezing line
            ax.plot(PSbins, PTFreezing, linestyle='--', linewidth=1., color='g', label='Freezing Line', zorder=4)

            # --- CHANGED: Replace Scatter Plot with a Volume-Weighted 2D Histogram ---
            # 150x150 bins provides an optimal trade-off between grid density and visual resolution
            counts, xedges, yedges, im = ax.hist2d(
                PS_flat, PT_flat,
                bins=[150, 150],
                range=[x_lim, y_lim],
                weights=Vol_flat,
                cmap='cmo.deep',
                norm=LogNorm(),
                cmin=1e-10,  # Do not draw bins containing zero volume
                zorder=1,
                rasterized=True
            )

            # Append a colorbar tracking total cell volumes per TS bin
            cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.04)
            cbar.set_label('Grid Cell Volume ($m^3$)', fontsize=11)
            # -------------------------------------------------------------------------

            # Core volume integrated centroid marker
            ax.plot(PS_core_avg, PT_core_avg, color='maroon', linestyle='None', marker='s', markersize=6, mec='k',
                    label='Vol-Weighted Mean', zorder=5)

            ax.set_ylim(y_lim)
            ax.set_xlim(x_lim)
            ax.set_xlabel('Salinity (PSU)', fontsize=12)
            ax.set_ylabel('Potential Temperature ($^\circ$C)', fontsize=12)
            ax.set_title(f"{region_name} | {dx}_{combined_cases_str}\n{years_str} ({label})", fontsize=11)
            ax.legend(loc='upper left', fontsize=8)

            plt.tight_layout()

            # Output filenames dynamically tracking the active averaging scheme
            out_filename = f"{dir_fig_save}/TS_{region_name}_{dx}_{combined_cases_str}_{years_str}_{label}.png"
            plt.savefig(out_filename, bbox_inches='tight', dpi=400)

            # Explicit figure flushing to avoid cross-process leaks
            fig.clear()
            plt.close(fig)
            print(f"--> [{dx}_{combined_cases_str}][{label}] Saved TS image to: {out_filename}")


# =========================================================================
# 2. Main Execution Orchestrator
# =========================================================================
if __name__ == "__main__":

    # -----------------------------------------------------------------
    # SPECIFY CONFIGURATIONS (From plot_spatial_stats.py)
    # -----------------------------------------------------------------
    RUN_TYPE = 'test' # e.g., 'Spin1' 'Spin6'
    TARGET_YEARS = ['0002', '0003', '0004']  # e.g., ['0002', '0003', '0004']
    regions_to_plot = ["FRIS", "RonneDcavity", "FilchnerDcavity", "RonneDshelf", "FilchnerDshelf", "BerknerBank",
                       "BerknerSouth", "FRISshelf"]  # Keys matching gmask_reg.py

    # Define the requested seasonal intervals.
    # Use lists of month integers. Set value to None to process all 12 months as Annual.
    seasonal_windows = {
        "Annual": None,
        "JFM": [1, 2, 3],
        "AMJ": [4, 5, 6],
        "JAS": [7, 8, 9],
        "OND": [11, 10, 12]
    }

    if RUN_TYPE == 'Spin1':
        simulations = [
            ('8', [('Spin1', 'p1')]),
            ('4', [('Spin1', 'p1')]),
            ('2', [('Spin1', 'p1'), ('Spin1', 'p2')]),
            ('1', [('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')])
        ]
    elif RUN_TYPE == 'Spin6':
        simulations = [
            ('8', [('Spin6', 'p1')]),
            ('8', [('Spin6', 'GMF1')]),
            ('4', [('Spin6', 'p1')]),
            ('2', [('Spin6', 'p1')]),
            ('1', [('Spin6', 'p1')])
        ]

    if RUN_TYPE == 'test':
        seasonal_windows = {
            "JFM": [1, 2, 3]
        }
        simulations = [
            ('8', [('Spin6', 'p1')])
        ]

    dir_fig_save = '/pscratch/sd/v/vankova/fris_analysis/fris_plots/TS_diagrams_histo'
    os.makedirs(dir_fig_save, exist_ok=True)

    # -----------------------------------------------------------------
    # PRE-COMPUTE SHARED BACKGROUND BACKGROUND DENSITY MATRIX VARIABLES
    # -----------------------------------------------------------------
    y_lim = np.array([-3.0, 1.0])
    x_lim = np.array([34.0, 35.2])
    PTbins = np.linspace(-3.5, 4, num=200)
    PSbins = np.linspace(32.0, 35.5, num=200)
    SAbins = gsw.SA_from_SP(PSbins, p=0., lon=0., lat=-75.)
    CTbins = gsw.pt_from_CT(SAbins, PTbins)
    CTgrid, SAgrid = np.meshgrid(CTbins, SAbins)
    PSgrid = gsw.SP_from_SA(SAgrid, p=0., lon=0., lat=-75.)
    PTgrid = gsw.pt_from_CT(SAgrid, CTgrid)
    neutralDensity = gsw.sigma0(SAgrid, CTgrid)
    rhoInterval = 0.2
    contours = np.arange(23., 29. + rhoInterval, rhoInterval)
    CTFreezing = gsw.CT_freezing(SAbins, 0, 1)
    PTFreezing = gsw.t_from_CT(SAbins, CTFreezing, p=0.)

    TS_bg_config = {
        'y_lim': y_lim, 'x_lim': x_lim, 'PSbins': PSbins,
        'PSgrid': PSgrid, 'PTgrid': PTgrid, 'neutralDensity': neutralDensity,
        'contours': contours, 'PTFreezing': PTFreezing
    }

    # -----------------------------------------------------------------
    # GENERATE PARALLEL TASKS BUNDLES
    # -----------------------------------------------------------------
    tasks = []
    for Fnum, cases in simulations:
        tasks.append(
            (Fnum, cases, RUN_TYPE, TARGET_YEARS, regions_to_plot, dir_fig_save, TS_bg_config, seasonal_windows)
        )

    # Allocate process thread pool based on work list volume (Up to 16 combinations)
    num_processes = min(len(tasks), 16)

    print(f"Spawning an isolated execution pool of {num_processes} parallel processes to generate TS diagrams...")

    with Pool(processes=num_processes) as pool:
        pool.map(process_single_ts_task, tasks)

    print("All parallel simulations and regional seasonal TS diagrams completed successfully.")