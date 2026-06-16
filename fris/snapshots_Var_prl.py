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
# 1. Worker Function for Parallel Processing
# =========================================================================
def process_single_file(file_path, idx, total_files, dx, sec, subsec_str, max_level_cell, opt_region, iis, lon, lat,
                        out_plot_dir, dsMesh_trimmed, var_config):
    """Processes a single history file and generates a plot based on variable configuration."""
    base_filename = os.path.basename(file_path)

    try:
        date_str = base_filename.split('.')[-2]
    except IndexError:
        date_str = f"idx_{idx:03d}"

    # Print with process ID to track progress in logs
    print(f"[{os.getpid()}] Processing {idx + 1}/{total_files}: {base_filename}")

    var_name = var_config['name']
    opt_proj = var_config['opt_proj']

    with xr.open_dataset(file_path) as ds:
        # Check if calculating the custom composite GMkappa variable
        if var_name == 'GMkappa':
            required_vars = ['timeMonthly_avg_gmBolusKappa', 'timeMonthly_avg_gmHorizontalTaper',
                             'timeMonthly_avg_gmKappaScaling']
            if all(v in ds for v in required_vars):
                # Pull top layer (index 0) for gmKappaScaling due to vertical scaling requirement
                scaling_top = ds['timeMonthly_avg_gmKappaScaling'].isel(Time=0, nVertLevelsP1=0).values
                bolus = ds['timeMonthly_avg_gmBolusKappa'].isel(Time=0).values
                taper = ds['timeMonthly_avg_gmHorizontalTaper'].isel(Time=0).values

                data_raw = (bolus * taper * scaling_top) / 1800.0
            else:
                missing = [v for v in required_vars if v not in ds]
                return f"Error: Missing variables {missing} for GMkappa calculation in {base_filename}"
        # Check if requested standard variable exists in dataset
        elif var_name in ds:
            data_raw = ds[var_name].isel(Time=0).values
        else:
            return f"Error: Variable '{var_name}' missing in {base_filename}"

    # Determine dimensions to selectively handle 3D vs 2D arrays
    if data_raw.ndim == 2:  # Shape: (nVertLevels, nCells)
        # If it has 2 dimensions, extract the bottom level using max_level_cell
        num_cells = data_raw.shape[1] if data_raw.shape[0] != max_level_cell.shape[0] else data_raw.shape[0]
        # Assuming shape is (nCells, nVertLevels) based on typical MPAS structure
        if data_raw.shape[0] == max_level_cell.shape[0]:
            plot_data = data_raw[np.arange(num_cells), max_level_cell]
        else:
            plot_data = data_raw[max_level_cell, np.arange(num_cells)]
        plot_data = np.where(max_level_cell >= 0, plot_data, np.nan)
    else:  # Shape: (nCells,)
        # Already 2D spatial field (e.g., MLD or Column Integrated Speed)
        plot_data = data_raw

    if opt_region:
        plot_data = np.where(iis == 0, np.nan, plot_data)

    # Plot generation
    if opt_proj == 'll':
        projection = ccrs.PlateCarree()
    elif opt_proj == 'sps':
        projection = ccrs.SouthPolarStereo(central_longitude=-75)

    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True, subplot_kw={"projection": projection})

    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)

    # Dynamically apply manual vmin, vmax, and colormap configurations
    collection = mosaic.polypcolor(
        ax, descriptor, plot_data,
        norm=colors.Normalize(vmin=var_config['vmin'], vmax=var_config['vmax']),
        cmap=var_config['cmap']
    )

    # Overlay Contour Lines dynamically if specified
    valid_mask = ~np.isnan(plot_data)
    if var_config['contours'] and np.any(valid_mask):
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

    fig.colorbar(collection, fraction=0.1, shrink=0.5, label=var_config['cb_label'])
    ax.set_title(f"{var_config['title_prefix']} - {dx} {sec} {subsec_str} ({date_str})", fontsize=12, pad=10)

    # Output file matches variable file_prefix definition
    output_png_path = f'{out_plot_dir}/{var_config["file_prefix"]}_{dx}_{sec}{subsec_str}_{date_str}.png'
    plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
    plt.close(fig)

    # --- Added print statement to track where images are saved ---
    print(f"[{os.getpid()}] Saved image to: {output_png_path}")

    return None


# =========================================================================
# 2. Main Logic Flow
# =========================================================================
if __name__ == "__main__":
    # NERSC Perlmutter allocates CPUs via environment variables
    NUM_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", 32))

    fris_loc = '/pscratch/sd/v/vankova/fris_analysis/fris_plots'
    opt_region = False
    iceshelves = ["Shelf"]

    # -----------------------------------------------------------------
    # VARIABLE CONFIGURATION DICTIONARY
    # Switch out, edit configurations, or loop over multiple dictionaries here.
    # -----------------------------------------------------------------

    PLOT_VARIABLE = 'GMkappa'  # Options: 'Tbot', 'ColSpeed', 'MLD', or 'GMkappa'

    if PLOT_VARIABLE == 'Tbot':
        VAR_CONFIG = {
            'name': 'timeMonthly_avg_activeTracers_temperature',  # Variable name inside netCDF
            'vmin': -2.6,
            'vmax': -1.8,
            'contours': [-1.9],  # List of contour values (or None/[] to skip)
            'cmap': 'cmo.thermal',
            'cb_label': 'Sea Floor Temperature [°C]',
            'title_prefix': 'Bottom Temperature & -1.9°C Contour',
            'file_prefix': 'Tbot' , # Used for naming the output file and directories
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
            'vmax': 1.0,  # Adjust limits based on your expected normalized values
            'contours': [],
            'cmap': 'CMRmap',
            'cb_label': 'GM Kappa / max(GM Kappa)',
            'title_prefix': 'GM Kappa',
            'file_prefix': 'GMkappa',
            'opt_proj': 'sps'
        }

    simulations = {
        '8': [('Spin6', 'p1')],
        '4': [('Spin6', 'p1')],
        '2': [('Spin6', 'p1')],
        '1': [('Spin6', 'p1')]
    }

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

            # Dynamically use the variable file prefix for the destination output folder
            out_plot_dir = f'{fris_loc}/snapshots_{VAR_CONFIG["file_prefix"]}/{dx}_{sec}{subsec_str}'
            os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

            file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"
            file_list = sorted(glob.glob(file_pattern))
            total_files = len(file_list)

            if not file_list:
                print(f"--> Warning: Missing files matching: {file_pattern}. Skipping.")
                continue

            print(
                f"\nSpinning up parallel pool with {NUM_WORKERS} workers for {total_files} files plotting '{VAR_CONFIG['name']}'...")

            # Pass the configuration dictionary right into the worker task
            worker_task = partial(
                process_single_file,
                total_files=total_files, dx=dx, sec=sec, subsec_str=subsec_str,
                max_level_cell=max_level_cell, opt_region=opt_region, iis=iis,
                lon=lon, lat=lat, out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed,
                var_config=VAR_CONFIG
            )

            iterable_args = [(file_path, idx) for idx, file_path in enumerate(file_list)]

            with mp.Pool(processes=NUM_WORKERS) as pool:
                results = pool.starmap(worker_task, iterable_args)

            for res in results:
                if res is not None:
                    print(res)

            print(f"Completed rendering all frames for configuration: {dx}_{sec}{subsec_str}.")