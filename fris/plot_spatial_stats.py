#!/usr/bin/env python3

import glob
import os
import numpy as np
import xarray as xr
import multiprocessing as mp
from functools import partial

import matplotlib
matplotlib.use('Agg')  # Force non-interactive backend (crucial for clusters)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean
import mosaic
from matplotlib import colors


# =========================================================================
# 1. Plotting Function for Statistical Fields
# =========================================================================
def generate_spatial_plot(plot_data, date_str, stat_type, dx, sec, subsec_str, lon, lat,
                          out_plot_dir, dsMesh_trimmed, var_config):
    """
    Generates a plot for a computed statistical field.
    stat_type: 'Avg', 'StdDev', 'MAD', or 'StdDevDiff'
    """
    opt_proj = var_config['opt_proj']

    if opt_proj == 'll':
        projection = ccrs.PlateCarree()
    elif opt_proj == 'sps':
        projection = ccrs.SouthPolarStereo(central_longitude=-75)

    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True, subplot_kw={"projection": projection})
    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)

    # Dynamic styling configurations based on statistical metric
    if stat_type in ['StdDev', 'MAD', 'StdDevDiff']:
        # Variability/dispersion fields are strictly positive.
        # Scale colorbar from 0 up to 50% of the standard variable range to capture fine fluctuations.
        vmin, vmax = 0.0, var_config['vmax'] * 0.5
        cmap = 'cmo.amp'
        draw_contours = False

        if stat_type == 'StdDev':
            cb_label = f"Std Dev: {var_config['cb_label']}"
            title_suffix = f"Standard Deviation ({date_str})"
            file_suffix = f"StdDev_{date_str}"
        elif stat_type == 'MAD':
            cb_label = f"MAD: {var_config['cb_label']}"
            title_suffix = f"Median Absolute Deviation ({date_str})"
            file_suffix = f"MAD_{date_str}"
        else:  # StdDevDiff
            cb_label = f"First-Diff Std Dev: {var_config['cb_label']}"
            title_suffix = f"Std Dev of Differenced Data ({date_str})"
            file_suffix = f"StdDevDiff_{date_str}"
    else:
        # Standard average field settings
        vmin, vmax = var_config['vmin'], var_config['vmax']
        cmap = var_config['cmap']
        cb_label = var_config['cb_label']
        title_suffix = f"Avg ({date_str})"
        file_suffix = f"Avg_{date_str}"
        draw_contours = True

    # Render spatial grid cells
    collection = mosaic.polypcolor(
        ax, descriptor, plot_data,
        norm=colors.Normalize(vmin=vmin, vmax=vmax),
        cmap=cmap,
        edgecolors='face',
        antialiased=False
    )

    # Overlay Contour Lines if plotting Mean
    valid_mask = ~np.isnan(plot_data)
    if draw_contours and var_config['contours'] and np.any(valid_mask):
        ax.tricontour(
            lon[valid_mask], lat[valid_mask], plot_data[valid_mask],
            levels=var_config['contours'], colors='black', linewidths=1.5, transform=ccrs.PlateCarree()
        )

    if opt_proj == 'll':
        ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())
    elif opt_proj == 'sps':
        ax.set_extent([-80, 0, -84, -64], ccrs.PlateCarree())

    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

    fig.colorbar(collection, fraction=0.1, shrink=1.0, label=cb_label)
    ax.set_title(f"{var_config['title_prefix']} - {title_suffix}\n{dx} {sec} {subsec_str}", fontsize=12, pad=10)

    # Output path construction
    output_png_path = f'{out_plot_dir}/{var_config["file_prefix"]}_{dx}_{sec}{subsec_str}_{file_suffix}.png'
    plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
    plt.close(fig)

    print(f"[{stat_type}] Saved image to: {output_png_path}")


# =========================================================================
# 2. Main Logic Flow
# =========================================================================
if __name__ == "__main__":
    NUM_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 32))

    fris_loc = '/pscratch/sd/v/vankova/fris_analysis/fris_plots/spatial_stats'
    opt_region = False
    iceshelves = ["Shelf"]

    # -----------------------------------------------------------------
    # SPECIFY TARGET YEARS FOR STATISTICS
    # -----------------------------------------------------------------
    RUN_TYPE = 'Spin6'
    TARGET_YEARS = ['0002', '0003', '0004']  # Example: Spin6 Years 2-4

    # Options: 'Tbot', 'ColSpeed', 'MLD', 'GMkappa', 'Melt', 'Ustar', or 'BLTemp'
    PLOT_VARIABLE = 'Melt'

    if RUN_TYPE == 'Spin1':
        simulations = {
            '8': [('Spin1', 'p1')],
            '4': [('Spin1', 'p1')],
            '2': [('Spin1', 'p1'), ('Spin1', 'p2')],
            '1': [('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
        }
    elif RUN_TYPE == 'Spin6':
        simulations = {
            '8': [('Spin6', 'p1')],
            '4': [('Spin6', 'p1')],
            '2': [('Spin6', 'p1')],
            '1': [('Spin6', 'p1')]
        }

    # -----------------------------------------------------------------
    # VARIABLE CONFIGURATIONS
    # -----------------------------------------------------------------
    if PLOT_VARIABLE == 'Tbot':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_activeTracers_temperature',
            'vmin': -2.6,
            'vmax': -1.8,
            'contours': [-1.9],
            'cmap': 'cmo.thermal',
            'cb_label': 'Sea Floor Temperature [°C]',
            'title_prefix': 'Bottom Temperature & -1.9°C Contour',
            'file_prefix': 'Tbot',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'ColSpeed':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_columnIntegratedSpeed',
            'vmin': 0.0,
            'vmax': 0.3,
            'contours': [],
            'cmap': 'cmo.speed',
            'cb_label': 'Column Integrated Speed [m/s]',
            'title_prefix': 'Column Integrated Speed',
            'file_prefix': 'ColSpeed',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'MLD':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_dThreshMLD',
            'vmin': 0,
            'vmax': 500,
            'contours': [],
            'cmap': 'cmo.deep',
            'cb_label': 'Mixed Layer Depth [m]',
            'title_prefix': 'MLD',
            'file_prefix': 'MLD',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'GMkappa':
        VAR_CONFIG = {
            'name': 'GMkappa',
            'vmin': 0.0,
            'vmax': 1.0,
            'contours': [],
            'cmap': 'CMRmap',
            'cb_label': 'GM Kappa / max(GM Kappa)',
            'title_prefix': 'GM Kappa',
            'file_prefix': 'GMkappa',
            'opt_proj': 'sps'
        }
    elif PLOT_VARIABLE == 'Melt':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_landIceFreshwaterFlux',
            'vmin': -3.0,
            'vmax': 3.0,      # Adjust based on expected maximum values
            'contours': [],
            'cmap': 'cmo.inflection',
            'cb_label': 'Melt rate interfacial [m/a]',
            'title_prefix': 'Melt rate Interfacial',
            'file_prefix': 'Melt',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'MeltTotal':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_landIceFreshwaterFluxTotal',
            'vmin': -3.0,
            'vmax': 3.0,      # Adjust based on expected maximum values
            'contours': [],
            'cmap': 'cmo.inflection',
            'cb_label': 'Melt rate total [m/a]',
            'title_prefix': 'Melt rate total',
            'file_prefix': 'MeltTot',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'Ustar':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_landIceFrictionVelocity',
            'vmin': 0.0,
            'vmax': 10.0,     # Adjust based on expected maximum values
            'contours': [],
            'cmap': 'cmo.speed',
            'cb_label': 'Ustar [cm/s]',
            'title_prefix': 'Land Ice Friction Velocity (Ustar)',
            'file_prefix': 'Ustar',
            'opt_proj': 'll'
        }
    elif PLOT_VARIABLE == 'BLTemp':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature',
            'vmin': -2.5,
            'vmax': 1.0,      # Adjust based on expected temperature bounds
            'contours': [],
            'cmap': 'cmo.thermal',
            'cb_label': 'BL temperature [°C]',
            'title_prefix': 'Boundary Layer Temperature',
            'file_prefix': 'BLTemp',
            'opt_proj': 'll'
        }

    # -----------------------------------------------------------------
    # CORE PROCESSING LOOP
    # -----------------------------------------------------------------
    for Fnum, cases in simulations.items():
        dx = f'F{Fnum}'
        run_name_mask = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
        fpath_mask = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name_mask}/run'
        mesh_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

        if not os.path.exists(mesh_file):
            print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}.")
            continue

        with xr.open_dataset(mesh_file) as dsMesh:
            dsMesh_trimmed = dsMesh[['latCell', 'lonCell', 'landIceFloatingMask', 'areaCell',
                                     'maxLevelCell', 'nEdges', 'nVertices', 'lonEdge', 'latEdge',
                                     'lonVertex', 'latVertex', 'cellsOnEdge', 'cellsOnVertex',
                                     'verticesOnEdge', 'verticesOnCell', 'edgesOnVertex']].load()

        lat = dsMesh_trimmed['latCell'].values * 180 / np.pi
        lon = dsMesh_trimmed['lonCell'].values * 180 / np.pi
        max_level_cell = dsMesh_trimmed['maxLevelCell'].values - 1

        iis = None
        if opt_region:
            import gmask_reg
            iam = gmask_reg.get_mask(iceshelves, mesh_file, opt_noGL=0, opt_wct=1)
            iis = iam[0, :]

        for sec, subsec in cases:
            subsec_str = subsec if sec == 'Spin1' else ''

            if sec == 'Spin6':
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

            out_plot_dir = f'{fris_loc}/statistics_{VAR_CONFIG["file_prefix"]}/{dx}_{sec}{subsec_str}'
            os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

            file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"
            all_files = sorted(glob.glob(file_pattern))

            # --- Filter target years and deduplicate monthly tracks ---
            unique_months_dict = {}
            for fp in all_files:
                try:
                    date_part = os.path.basename(fp).split('.')[-2]
                    ymd = date_part.split('_')[0]
                    file_year, file_month, _ = ymd.split('-')

                    if file_year in TARGET_YEARS:
                        unique_months_dict[(file_year, file_month)] = fp
                except (IndexError, ValueError):
                    continue

            filtered_files = [unique_months_dict[k] for k in sorted(unique_months_dict.keys())]

            # --- Strict File Count Verification Check ---
            expected_count = 12 * len(TARGET_YEARS)
            actual_count = len(filtered_files)

            if actual_count != expected_count:
                raise ValueError(
                    f"\n[ERROR] File count mismatch for {dx}_{sec}{subsec_str}!\n"
                    f"Expected exactly {expected_count} files ({len(TARGET_YEARS)} years * 12 months), "
                    f"but found {actual_count} files.\n"
                    f"Target Years: {TARGET_YEARS}\n"
                )

            print(f"\nProcessing time statistics over years {TARGET_YEARS} using exactly {len(filtered_files)} files...")

            var_name = VAR_CONFIG['name']

            # --- Open and Process Timeseries Statistics ---
            with xr.open_mfdataset(filtered_files, chunks={'Time': 1}, combine='nested', concat_dim='Time') as ds:
                if var_name == 'GMkappa':
                    required_vars = ['timeMonthly_avg_gmBolusKappa', 'timeMonthly_avg_gmHorizontalTaper',
                                     'timeMonthly_avg_gmKappaScaling']
                    if all(v in ds for v in required_vars):
                        scaling_top = ds['timeMonthly_avg_gmKappaScaling'].isel(nVertLevelsP1=0)
                        bolus = ds['timeMonthly_avg_gmBolusKappa']
                        taper = ds['timeMonthly_avg_gmHorizontalTaper']

                        ts_data = (bolus * taper * scaling_top) / 1800.0
                    else:
                        print(f"Error: Missing variables for GMkappa calculations.")
                        continue
                elif var_name in ds:
                    ts_data = ds[var_name]
                else:
                    print(f"Error: Variable '{var_name}' missing in dataset.")
                    continue

                # --- Apply conversions based on the current field parameter ---
                if var_name == 'timeMonthly_avg_landIceFreshwaterFlux' or var_name == 'timeMonthly_avg_landIceFreshwaterFluxTotal':
                    sec_per_year = 365.25 * 24 * 3600  # ~31,557,600 seconds
                    rho_fw = 1000.0                   # Freshwater density (kg/m^3)
                    ts_data = ts_data * sec_per_year / rho_fw

                elif var_name == 'timeMonthly_avg_landIceFrictionVelocity':
                    m2cm = 100
                    ts_data = ts_data * m2cm

                # Compute core statistical operations over the 'Time' dimension
                raw_mean = ts_data.mean(dim='Time').values
                raw_std = ts_data.std(dim='Time').values

                # 1. Median Absolute Deviation calculation: median(|x - median(x)|)
                ts_median = ts_data.median(dim='Time')
                raw_mad = np.abs(ts_data - ts_median).median(dim='Time').values

                # 2. Standard Deviation of Differenced Data (diff over Time axis)
                raw_std_diff = ts_data.diff(dim='Time').std(dim='Time').values

            # Process 3D structural slices down to bottom 2D topologies across statistics
            stats_to_plot = [
                ('Avg', raw_mean),
                ('StdDev', raw_std),
                ('MAD', raw_mad),
                ('StdDevDiff', raw_std_diff)
            ]
            processed_fields = {}

            for label, data_raw in stats_to_plot:
                if data_raw.ndim == 2:
                    num_cells = data_raw.shape[1] if data_raw.shape[0] != max_level_cell.shape[0] else data_raw.shape[0]
                    if data_raw.shape[0] == max_level_cell.shape[0]:
                        plot_data = data_raw[np.arange(num_cells), max_level_cell]
                    else:
                        plot_data = data_raw[max_level_cell, np.arange(num_cells)]
                    plot_data = np.where(max_level_cell >= 0, plot_data, np.nan)
                else:
                    plot_data = data_raw

                if opt_region:
                    plot_data = np.where(iis == 0, np.nan, plot_data)

                processed_fields[label] = plot_data

            # Define historical tag markers for output filenames
            if len(TARGET_YEARS) == 1:
                date_tag = f"year_{TARGET_YEARS[0]}"
            else:
                date_tag = f"years_{TARGET_YEARS[0]}-{TARGET_YEARS[-1]}"

            # --- Generate the 4 Spatial Plots ---
            for metric in ['Avg', 'StdDev', 'MAD', 'StdDevDiff']:
                generate_spatial_plot(
                    plot_data=processed_fields[metric], date_str=date_tag, stat_type=metric,
                    dx=dx, sec=sec, subsec_str=subsec_str, lon=lon, lat=lat,
                    out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed, var_config=VAR_CONFIG
                )

            print(f"Completed rendering all statistical fields for: {dx}_{sec}{subsec_str}.\n")