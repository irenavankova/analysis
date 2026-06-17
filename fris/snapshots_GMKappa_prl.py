#!/usr/bin/env python3

import glob
import os
import numpy as np
import xarray as xr

import matplotlib

matplotlib.use('Agg')  # Force non-interactive backend (crucial for clusters)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mosaic
from matplotlib import colors


# =========================================================================
# 1. Processing and Plotting Function
# =========================================================================
def generate_annual_plot(plot_data, target_year, dx, cases_str, max_level_cell, opt_region, iis, lon, lat,
                         out_plot_dir, dsMesh_trimmed, var_config):
    """Generates an annual average plot based on averaged raw cell data."""

    opt_proj = var_config['opt_proj']

    # Determine dimensions to selectively handle 3D vs 2D arrays
    if plot_data.ndim == 2:  # Shape: (nVertLevels, nCells) or (nCells, nVertLevels)
        num_cells = plot_data.shape[1] if plot_data.shape[0] != max_level_cell.shape[0] else plot_data.shape[0]
        if plot_data.shape[0] == max_level_cell.shape[0]:
            spatial_data = plot_data[np.arange(num_cells), max_level_cell]
        else:
            spatial_data = plot_data[max_level_cell, np.arange(num_cells)]
        spatial_data = np.where(max_level_cell >= 0, spatial_data, np.nan)
    else:  # Shape: (nCells,)
        spatial_data = plot_data

    if opt_region and iis is not None:
        spatial_data = np.where(iis == 0, np.nan, spatial_data)

    # Plot generation
    if opt_proj == 'll':
        projection = ccrs.PlateCarree()
    elif opt_proj == 'sps':
        projection = ccrs.SouthPolarStereo(central_longitude=-75)
    else:
        projection = ccrs.PlateCarree()

    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True, subplot_kw={"projection": projection})

    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)

    # Dynamically apply manual vmin, vmax, and colormap configurations
    collection = mosaic.polypcolor(
        ax, descriptor, spatial_data,
        norm=colors.Normalize(vmin=var_config['vmin'], vmax=var_config['vmax']),
        cmap=var_config['cmap']
    )

    # Overlay Contour Lines dynamically if specified
    valid_mask = ~np.isnan(spatial_data)
    if var_config['contours'] and np.any(valid_mask):
        ax.tricontour(
            lon[valid_mask], lat[valid_mask], spatial_data[valid_mask],
            levels=var_config['contours'], colors='black', linewidths=1.5, transform=ccrs.PlateCarree()
        )

    if opt_proj == 'll':
        ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())
    elif opt_proj == 'sps':
        ax.set_extent([-80, 0, -84, -64], ccrs.PlateCarree())

    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

    fig.colorbar(collection, fraction=0.1, shrink=1.0, label=var_config['cb_label'])
    ax.set_title(f"{var_config['title_prefix']} (Annual Avg Year {target_year:04d}) - {dx} {cases_str}", fontsize=12,
                 pad=10)

    # Output file configuration
    date_str = f"Year{target_year:04d}_AnnualAvg"
    output_png_path = f'{out_plot_dir}/{var_config["file_prefix"]}_{dx}_{cases_str}_{date_str}.png'
    plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
    plt.close(fig)

    print(f"--> Saved annual average image to: {output_png_path}\n")


# =========================================================================
# 2. Main Logic Flow
# =========================================================================
if __name__ == "__main__":

    fris_loc = '/pscratch/sd/v/vankova/fris_analysis/fris_plots'
    opt_region = False
    iceshelves = ["Shelf"]

    # --- TOGGLE VARIABLES ---
    TARGET_YEAR = 5  # Set the desired simulation year to aggregate and average

    VAR_CONFIG = {
        'name': 'GMkappa',
        'vmin': 0.0,
        'vmax': 1.0,
        'contours': [],
        'cmap': 'CMRmap_r',
        'cb_label': 'GM Kappa / max(GM Kappa)',
        'title_prefix': 'GM Kappa Annual Average',
        'file_prefix': 'GMkappa_Annual',
        'opt_proj': 'sps'
    }

    simulations = {
        '8': [('Spin1', 'p1')],
        '4': [('Spin1', 'p1')],
        '2': [('Spin1', 'p1'), ('Spin1', 'p2')],
        '1': [('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
    }

    for Fnum, cases in simulations.items():
        dx = f'F{Fnum}'
        print(f"=========================================================================")
        print(f" Processing Resolution: {dx}")
        print(f"=========================================================================")

        # Use the first available setup to grab geometry mesh data
        run_name_mask = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
        fpath_mask = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name_mask}/run'
        mesh_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

        if not os.path.exists(mesh_file):
            print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}. Skipping resolution.")
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

        # Gather files for Year 5 across ALL folders/sub-sequences for this resolution
        year_file_list = []
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

            # Search strictly for files matching the requested target year (e.g., .0005-*-*.nc)
            file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.{TARGET_YEAR:04d}-*-*.nc"
            found_files = glob.glob(file_pattern)
            year_file_list.extend(found_files)

        # Sort files to ensure chronological processing
        year_file_list = sorted(list(set(year_file_list)))

        if not year_file_list:
            print(f"--> Warning: No files found matching Year {TARGET_YEAR:04d} for {dx}. Skipping.")
            continue

        print(f"Files selected for Year {TARGET_YEAR:04d} Annual Average ({dx}):")
        for f in year_file_list:
            print(f"  - {os.path.basename(f)}")
        print(f"Total files found: {len(year_file_list)}")

        # Accumulator array list for calculation
        monthly_gm_kappa_arrays = []

        # Process gathered netCDF datasets sequentially to calculate and average fields
        for file_path in year_file_list:
            with xr.open_dataset(file_path) as ds:
                required_vars = ['timeMonthly_avg_gmBolusKappa', 'timeMonthly_avg_gmHorizontalTaper',
                                 'timeMonthly_avg_gmKappaScaling']

                if all(v in ds for v in required_vars):
                    # Compute GMkappa spatial field slice
                    scaling_top = ds['timeMonthly_avg_gmKappaScaling'].isel(Time=0, nVertLevelsP1=0).values
                    bolus = ds['timeMonthly_avg_gmBolusKappa'].isel(Time=0).values
                    taper = ds['timeMonthly_avg_gmHorizontalTaper'].isel(Time=0).values

                    gm_kappa_month = (bolus * taper * scaling_top) / 1800.0
                    monthly_gm_kappa_arrays.append(gm_kappa_month)
                else:
                    missing = [v for v in required_vars if v not in ds]
                    print(f"--> Error: Missing variables {missing} in {os.path.basename(file_path)}. Skipping file.")

        # --- STRICT 12-MONTH VALIDATION CHECK ---
        num_valid_months = len(monthly_gm_kappa_arrays)
        if num_valid_months != 12:
            raise ValueError(
                f"CRITICAL ERROR: Annual average calculation failed for resolution {dx}. "
                f"Expected exactly 12 valid monthly files for Year {TARGET_YEAR:04d}, "
                f"but found and processed {num_valid_months} files."
            )

        # Perform the actual annual mean average calculations across months
        annual_avg_gm_kappa = np.mean(np.array(monthly_gm_kappa_arrays), axis=0)

        # Combine strings of cases included for directory/plot titles (e.g. Spin1p1_Spin1p2)
        combined_cases_str = "_".join(cases_processed)
        out_plot_dir = f'{fris_loc}/snapshots_{VAR_CONFIG["file_prefix"]}_annual/{dx}_{combined_cases_str}'
        os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

        # Plot the consolidated map structure
        generate_annual_plot(
            plot_data=annual_avg_gm_kappa, target_year=TARGET_YEAR, dx=dx,
            cases_str=combined_cases_str, max_level_cell=max_level_cell,
            opt_region=opt_region, iis=iis, lon=lon, lat=lat,
            out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed,
            var_config=VAR_CONFIG
        )