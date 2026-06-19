#!/usr/bin/env python3

import glob
import os
import numpy as np
import xarray as xr
from multiprocessing import Pool

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
def generate_spatial_plot(plot_data, date_str, stat_type, dx, cases_str, opt_region, iis, lon, lat,
                          out_plot_dir, dsMesh_trimmed, var_config):
    """
    Generates a plot for a computed statistical field of sea ice variables (strictly 2D/surface fields).
    """
    opt_proj = var_config['opt_proj']

    if opt_proj == 'll':
        projection = ccrs.PlateCarree()
        fig_size = (5, 4)
    elif opt_proj in ['fris', 'sps', 'wed']:
        projection = ccrs.SouthPolarStereo(central_longitude=-75)
        fig_size = (6, 4.5)
    else:
        projection = ccrs.PlateCarree()

    # Data is strictly 2D/surface mapping to (nCells,)
    spatial_data = plot_data

    if opt_region and iis is not None:
        spatial_data = np.where(iis == 0, np.nan, spatial_data)

    fig, ax = plt.subplots(figsize=fig_size, constrained_layout=True, subplot_kw={"projection": projection})
    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)

    # Dynamic styling configurations based on statistical metric
    if stat_type in ['StdDev', 'MAD', 'StdDevDiff', 'StdDevWeighted']:
        if stat_type == 'StdDevWeighted':
            vmin, vmax = 0.0, 1.0
            cmap = 'CMRmap_r'
            cb_label = f"Std Dev / Mean (Dimensionless)"
            title_suffix = f"Mean-Weighted Std Dev ({date_str})"
            file_suffix = f"StdDevWeighted_{date_str}"
        else:
            range_span = abs(var_config['vmax'] - var_config['vmin'])
            scale_factor = var_config.get('vmax_scale_factor', 0.25)
            vmin, vmax = 0.0, range_span * scale_factor
            cmap = 'CMRmap_r'

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

        draw_contours = False
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
        ax, descriptor, spatial_data,
        norm=colors.Normalize(vmin=vmin, vmax=vmax),
        cmap=cmap,
        linewidth=0.0,
        edgecolors='face'
    )

    # Overlay Contour Lines if plotting Mean
    valid_mask = ~np.isnan(spatial_data)
    if draw_contours and var_config['contours'] and np.any(valid_mask):
        ax.tricontour(
            lon[valid_mask], lat[valid_mask], spatial_data[valid_mask],
            levels=var_config['contours'], colors='black', linewidths=1.5, transform=ccrs.PlateCarree()
        )

    if opt_proj == 'll':
        ax.set_extent([-85, -20, -84, -72], ccrs.PlateCarree())
    elif opt_proj == 'fris':
        ax.set_extent([-82, -25, -81, -74], ccrs.PlateCarree())
    elif opt_proj == 'wed':
        ax.set_extent([-82, -25, -81, -72], ccrs.PlateCarree())
    elif opt_proj == 'sps':
        ax.set_extent([-80, 0, -84, -64], ccrs.PlateCarree())

    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

    fig.colorbar(collection, fraction=0.1, shrink=1.0, label=cb_label, extend='both')

    # Output path construction
    output_png_path = f'{out_plot_dir}/{var_config["file_prefix"]}_{dx}_{cases_str}_{file_suffix}.png'
    plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
    plt.close(fig)

    print(f"--> [{stat_type}] Saved image to: {output_png_path}")


# =========================================================================
# 2. Worker Function for a Single Resolution
# =========================================================================
def process_single_resolution(args):
    """Processes a single spatial resolution key-value pair in a parallel worker."""
    Fnum, cases, VAR_CONFIG, RUN_TYPE, TARGET_YEARS, fris_loc, opt_region, iceshelves, PLOT_VARIABLE = args
    dx = f'F{Fnum}'

    print(f"=========================================================================\n"
          f" Starting Execution: {dx} | Variable: {PLOT_VARIABLE}\n"
          f"=========================================================================")

    run_name_mask = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
    fpath_mask = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name_mask}/run'
    mesh_file = f'{fpath_mask}/{run_name_mask}.mpassi.rst.0002-01-01_00000.nc'

    if not os.path.exists(mesh_file):
        print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}. Skipping resolution.")
        return

    with xr.open_dataset(mesh_file) as dsMesh:
        dsMesh_trimmed = dsMesh[['latCell', 'lonCell', 'areaCell',
                                 'nEdges', 'nVertices', 'lonEdge', 'latEdge',
                                 'lonVertex', 'latVertex', 'cellsOnEdge', 'cellsOnVertex',
                                 'verticesOnEdge', 'verticesOnCell', 'edgesOnVertex']].load()

    lat = dsMesh_trimmed['latCell'].values * 180 / np.pi
    lon = dsMesh_trimmed['lonCell'].values * 180 / np.pi

    iis = None
    if opt_region:
        import gmask_reg
        iam = gmask_reg.get_mask(iceshelves, mesh_file, opt_noGL=0, opt_wct=1)
        iis = iam[0, :]

    unique_months_dict = {}
    cases_processed = []

    for sec, subsec in cases:
        subsec_str = subsec if sec == 'Spin1' else ''
        cases_processed.append(f"{sec}{subsec_str}")

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

        for yr_str in TARGET_YEARS:
            try:
                yr_int = int(yr_str)
                file_pattern = f"{fpath}/{run_name}.mpassi.hist.am.timeSeriesStatsMonthly.{yr_int:04d}-*-*.nc"
                found_files = sorted(glob.glob(file_pattern))

                for file_path in found_files:
                    base_name = os.path.basename(file_path)
                    date_part = base_name.split('.')[-2]
                    unique_months_dict[date_part] = file_path
            except (IndexError, ValueError):
                continue

    year_file_list = [unique_months_dict[k] for k in sorted(unique_months_dict.keys())]

    if not year_file_list:
        print(f"--> Warning: No files found matching target years {TARGET_YEARS} for {dx}. Skipping.")
        return

    expected_count = 12 * len(TARGET_YEARS)
    num_valid_months = len(year_file_list)
    if num_valid_months != expected_count:
        raise ValueError(
            f"CRITICAL ERROR: Statistics calculation failed for resolution {dx}. "
            f"Expected exactly {expected_count} unique valid monthly files for Years {TARGET_YEARS}, "
            f"but found and processed {num_valid_months} unique files."
        )

    print(f"[{dx} | {PLOT_VARIABLE}] Files selected for Years {TARGET_YEARS} Processing:")
    for f in year_file_list:
        print(f"  - {os.path.basename(f)}")

    monthly_data_arrays = []
    skip_case = False

    # Scale factor to convert m/s to m/yr
    units_scale_factor = 60 * 60 * 24 * 365

    for file_path in year_file_list:
        with xr.open_dataset(file_path) as ds:
            try:
                if PLOT_VARIABLE == 'ice_concentration':
                    data_month = ds['timeMonthly_avg_iceAreaCell'].isel(Time=0).values
                elif PLOT_VARIABLE == 'ice_thickness':
                    data_month = ds['timeMonthly_avg_iceVolumeCell'].isel(Time=0).values
                elif PLOT_VARIABLE == 'ice_production':
                    congelation = ds['timeMonthly_avg_congelation'].isel(Time=0).values
                    frazil = ds['timeMonthly_avg_frazilFormation'].isel(Time=0).values
                    snowice = ds['timeMonthly_avg_snowiceFormation'].isel(Time=0).values
                    data_month = (congelation + frazil + snowice) * units_scale_factor
                elif PLOT_VARIABLE == 'ice_melting':
                    basal = ds['timeMonthly_avg_basalIceMelt'].isel(Time=0).values
                    surface = ds['timeMonthly_avg_surfaceIceMelt'].isel(Time=0).values
                    lateral = ds['timeMonthly_avg_lateralIceMelt'].isel(Time=0).values
                    data_month = (basal + surface + lateral) * units_scale_factor
                else:
                    print(f"--> Error: Requested PLOT_VARIABLE '{PLOT_VARIABLE}' extraction logic not set up.")
                    skip_case = True
                    break

                monthly_data_arrays.append(data_month)
            except KeyError as e:
                print(f"--> Error: Missing variable {e} in dataset {os.path.basename(file_path)}.")
                skip_case = True
                break

    if skip_case or not monthly_data_arrays:
        return

    data_all = np.array(monthly_data_arrays)

    raw_mean = np.mean(data_all, axis=0)
    raw_std = np.std(data_all, axis=0)

    data_median = np.median(data_all, axis=0)
    raw_mad = np.median(np.abs(data_all - data_median[np.newaxis, ...]), axis=0)

    data_diff = np.diff(data_all, axis=0)
    raw_std_diff = np.std(data_diff, axis=0)

    raw_std_weighted = np.where(raw_mean != 0, raw_std / raw_mean, np.nan)

    stats_to_plot = [
        ('Avg', raw_mean),
        ('StdDev', raw_std),
        ('MAD', raw_mad),
        ('StdDevDiff', raw_std_diff),
        ('StdDevWeighted', raw_std_weighted)
    ]

    combined_cases_str = "_".join(cases_processed)
    out_plot_dir = f'{fris_loc}/statistics_{VAR_CONFIG["file_prefix"]}/{dx}_{combined_cases_str}'
    os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

    if len(TARGET_YEARS) == 1:
        date_tag = f"Year{int(TARGET_YEARS[0]):04d}"
    else:
        date_tag = f"Years{int(TARGET_YEARS[0]):04d}-{int(TARGET_YEARS[-1]):04d}"

    for metric, field_data in stats_to_plot:
        generate_spatial_plot(
            plot_data=field_data, date_str=date_tag, stat_type=metric,
            dx=dx, cases_str=combined_cases_str,
            opt_region=opt_region, iis=iis, lon=lon, lat=lat,
            out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed, var_config=VAR_CONFIG
        )

    print(f"Finished processing and rendering all fields for: {dx}_{combined_cases_str} ({PLOT_VARIABLE}).\n")


# =========================================================================
# Helper function to dynamically retrieve configurations
# =========================================================================
def get_variable_config(variable_name):
    """Returns runtime mapping parameters optimized for sea ice configurations."""
    if variable_name == 'ice_concentration':
        return {
            'vmin': 0.0, 'vmax': 1.0, 'contours': [], 'cmap': 'cmo.ice',
            'cb_label': 'Ice Concentration [Fraction]', 'title_prefix': 'Sea Ice Concentration',
            'file_prefix': 'ice_concentration', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'ice_thickness':
        return {
            'vmin': 0.0, 'vmax': 4.0, 'contours': [], 'cmap': 'cmo.ice',
            'cb_label': 'Ice Thickness [m]', 'title_prefix': 'Sea Ice Thickness',
            'file_prefix': 'ice_thickness', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'ice_production':
        return {
            'vmin': -5.0, 'vmax': 5.0, 'contours': [], 'cmap': 'RdBu_r',
            'cb_label': 'Ice Production Rate [m/yr]', 'title_prefix': 'Total Sea Ice Production',
            'file_prefix': 'ice_production', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'ice_melting':
        return {
            'vmin': 0.0, 'vmax': 5.0, 'contours': [], 'cmap': 'cmo.amp',
            'cb_label': 'Ice Melting Rate [m/yr]', 'title_prefix': 'Total Sea Ice Melting',
            'file_prefix': 'ice_melting', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    else:
        raise ValueError(f"Unknown sea ice variable configuration requested: {variable_name}")


# =========================================================================
# 3. Main Execution Block
# =========================================================================
if __name__ == "__main__":

    fris_loc = '/pscratch/sd/v/vankova/fris_analysis/fris_plots/spatial_stats_si'
    opt_region = False
    iceshelves = ["Shelf"]

    # -----------------------------------------------------------------
    # SPECIFY TARGET SIMULATION TYPE AND YEARS
    # -----------------------------------------------------------------
    RUN_TYPE = 'Spin6'
    TARGET_YEARS = ['0002', '0003', '0004']

    # Array of target variables requested for Sea Ice analysis
    PLOT_VARIABLES = ['ice_concentration', 'ice_thickness', 'ice_production', 'ice_melting']

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
    # PARALLEL TASKS GENERATION
    # -----------------------------------------------------------------
    tasks = []
    for var in PLOT_VARIABLES:
        var_config = get_variable_config(var)
        for Fnum, cases in simulations.items():
            tasks.append(
                (Fnum, cases, var_config, RUN_TYPE, TARGET_YEARS, fris_loc, opt_region, iceshelves, var)
            )

    num_processes = min(len(tasks), 16)

    print(f"Spawning an isolated execution pool of {num_processes} parallel processes to process {len(tasks)} tasks...")

    with Pool(processes=num_processes) as pool:
        pool.map(process_single_resolution, tasks)

    print("All spatial resolution pools and multi-variable configurations completed successfully.")