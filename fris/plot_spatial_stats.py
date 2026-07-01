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
def generate_spatial_plot(plot_data, date_str, stat_type, dx, cases_str, max_level_cell, opt_region, iis, lon, lat,
                          out_plot_dir, dsMesh_trimmed, var_config):
    """
    Generates a plot for a computed statistical field.
    stat_type: 'Avg', 'StdDev', 'MAD', 'StdDevDiff', or 'StdDevWeighted'
    """
    opt_proj = var_config['opt_proj']

    if opt_proj == 'll':
        projection = ccrs.PlateCarree()
        fig_size = (5, 4)
    else:
        projection = ccrs.SouthPolarStereo(central_longitude=-75)
        fig_size = (6, 4.5)

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

    fig, ax = plt.subplots(figsize=fig_size, constrained_layout=True, subplot_kw={"projection": projection})
    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)

    # Dynamic styling configurations based on statistical metric
    if stat_type in ['StdDev', 'MAD', 'StdDevDiff', 'StdDevWeighted']:
        # Variability/dispersion fields are strictly positive.
        if stat_type == 'StdDevWeighted':
            # Relative variation metrics benefit from a fixed 0 to 1 scaling, or customized bounds
            vmin, vmax = 0.0, 1.0
            cmap = 'CMRmap_r'
            cb_label = f"Std Dev / Mean (Dimensionless)"
            title_suffix = f"Mean-Weighted Std Dev ({date_str})"
            file_suffix = f"StdDevWeighted_{date_str}"
        else:
            # Scale colorbar from 0 up to variable-specific range fraction to capture fine fluctuations.
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
        # antialiased=False
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
        # ax.set_extent([-85, -25, -84, -74], ccrs.PlateCarree())
        ax.set_extent([-82, -25, -81, -74], ccrs.PlateCarree())
    elif opt_proj == 'wed':
        ax.set_extent([-82, -25, -81, -72], ccrs.PlateCarree())
    elif opt_proj == 'sps':
        ax.set_extent([-80, 0, -84, -64], ccrs.PlateCarree())
    elif opt_proj == 'ross':
        ax.set_extent([165, 210, -86, -72], ccrs.PlateCarree())  # Ross with shelf

    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

    fig.colorbar(collection, fraction=0.1, shrink=1.0, label=cb_label, extend='both')
    # ax.set_title(f"{var_config['title_prefix']} - {title_suffix}\n{dx} {cases_str}", fontsize=12, pad=10)

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
    mesh_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

    if not os.path.exists(mesh_file):
        print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}. Skipping resolution.")
        return

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
                file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.{yr_int:04d}-*-*.nc"
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
    var_name = VAR_CONFIG['name']
    skip_case = False

    for file_path in year_file_list:
        with xr.open_dataset(file_path) as ds:
            if var_name in ['timeMonthly_avg_activeTracers_temperature', 'timeMonthly_avg_activeTracers_salinity']:
                if var_name in ds:
                    data_3d = ds[var_name].isel(Time=0).values
                    num_cells = data_3d.shape[0]

                    if PLOT_VARIABLE in ['Tint', 'Sint']:
                        if 'timeMonthly_avg_layerThickness' in ds:
                            h_3d = ds['timeMonthly_avg_layerThickness'].isel(Time=0).values

                            num_levels = data_3d.shape[1]
                            level_indices = np.arange(num_levels)[None, :]
                            valid_vertical_mask = level_indices <= max_level_cell[:, None]

                            h_masked = np.where(valid_vertical_mask, h_3d, 0.0)
                            data_masked = np.where(valid_vertical_mask, data_3d, 0.0)

                            total_thickness = np.sum(h_masked, axis=1)
                            weighted_sum = np.sum(data_masked * h_masked, axis=1)

                            averaged_data = np.where((total_thickness > 0) & (max_level_cell >= 0),
                                                     weighted_sum / total_thickness, np.nan)
                            monthly_data_arrays.append(averaged_data)
                        else:
                            print(
                                f"--> Error: 'timeMonthly_avg_layerThickness' missing for integration in {os.path.basename(file_path)}.")
                            skip_case = True
                            break
                    elif PLOT_VARIABLE in ['Tbot', 'Sbot']:
                        # Bottom layer selection
                        bottom_data = data_3d[np.arange(num_cells), max_level_cell]
                        bottom_data = np.where(max_level_cell >= 0, bottom_data, np.nan)
                        monthly_data_arrays.append(bottom_data)
                    elif PLOT_VARIABLE in ['Tsurf', 'Ssurf']:
                        # Top layer selection (index 0)
                        surf_data = data_3d[:, 0]
                        surf_data = np.where(max_level_cell >= 0, surf_data, np.nan)
                        monthly_data_arrays.append(surf_data)
                else:
                    print(f"--> Error: Variable '{var_name}' missing in dataset {os.path.basename(file_path)}.")
                    skip_case = True
                    break
            elif PLOT_VARIABLE in ['BotSpeed', 'SurfSpeed']:
                if 'timeMonthly_avg_velocityZonal' in ds and 'timeMonthly_avg_velocityMeridional' in ds:
                    u_3d = ds['timeMonthly_avg_velocityZonal'].isel(Time=0).values
                    v_3d = ds['timeMonthly_avg_velocityMeridional'].isel(Time=0).values

                    if PLOT_VARIABLE == 'BotSpeed':
                        num_cells = u_3d.shape[0]
                        u_target = u_3d[np.arange(num_cells), max_level_cell]
                        v_target = v_3d[np.arange(num_cells), max_level_cell]
                    else:  # SurfSpeed
                        u_target = u_3d[:, 0]
                        v_target = v_3d[:, 0]

                    speed = np.sqrt(u_target ** 2 + v_target ** 2)
                    speed = np.where(max_level_cell >= 0, speed, np.nan)
                    monthly_data_arrays.append(speed)
                else:
                    print(
                        f"--> Error: Zonal/Meridional velocity missing for {PLOT_VARIABLE} in {os.path.basename(file_path)}.")
                    skip_case = True
                    break
            elif PLOT_VARIABLE == 'DepthAvgSpeed':
                if all(v in ds for v in ['timeMonthly_avg_velocityZonal', 'timeMonthly_avg_velocityMeridional',
                                         'timeMonthly_avg_layerThickness']):
                    u_3d = ds['timeMonthly_avg_velocityZonal'].isel(Time=0).values
                    v_3d = ds['timeMonthly_avg_velocityMeridional'].isel(Time=0).values
                    h_3d = ds['timeMonthly_avg_layerThickness'].isel(Time=0).values

                    num_levels = u_3d.shape[1]
                    level_indices = np.arange(num_levels)[None, :]
                    valid_vertical_mask = level_indices <= max_level_cell[:, None]

                    h_masked = np.where(valid_vertical_mask, h_3d, 0.0)
                    total_thickness = np.sum(h_masked, axis=1)

                    u_masked = np.where(valid_vertical_mask, u_3d, 0.0)
                    v_masked = np.where(valid_vertical_mask, v_3d, 0.0)

                    u_avg = np.sum(u_masked * h_masked, axis=1) / np.where(total_thickness > 0, total_thickness, 1.0)
                    v_avg = np.sum(v_masked * h_masked, axis=1) / np.where(total_thickness > 0, total_thickness, 1.0)

                    speed_avg = np.sqrt(u_avg ** 2 + v_avg ** 2)
                    speed_avg = np.where((total_thickness > 0) & (max_level_cell >= 0), speed_avg, np.nan)
                    monthly_data_arrays.append(speed_avg)
                else:
                    print(
                        f"--> Error: Missing velocity or thickness keys for DepthAvgSpeed in {os.path.basename(file_path)}.")
                    skip_case = True
                    break
            elif var_name == 'Tstar':
                req_tstar = ['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature',
                             'timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']
                if all(v in ds for v in req_tstar):
                    tcb = ds['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature'].isel(
                        Time=0).values
                    tci = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature'].isel(Time=0).values
                    monthly_data_arrays.append(tcb - tci)
                else:
                    missing = [v for v in req_tstar if v not in ds]
                    print(f"--> Error: Missing variables {missing} for Tstar in {os.path.basename(file_path)}.")
                    skip_case = True
                    break
            elif var_name == 'GMkappa':
                required_vars = ['timeMonthly_avg_gmBolusKappa', 'timeMonthly_avg_gmHorizontalTaper',
                                 'timeMonthly_avg_gmKappaScaling']
                if all(v in ds for v in required_vars):
                    scaling_top = ds['timeMonthly_avg_gmKappaScaling'].isel(Time=0, nVertLevelsP1=0).values
                    bolus = ds['timeMonthly_avg_gmBolusKappa'].isel(Time=0).values
                    taper = ds['timeMonthly_avg_gmHorizontalTaper'].isel(Time=0).values

                    data_month = (bolus * taper * scaling_top) / 1800.0
                    monthly_data_arrays.append(data_month)
                else:
                    missing = [v for v in required_vars if v not in ds]
                    print(f"--> Error: Missing variables {missing} in {os.path.basename(file_path)}. Skipping.")
                    skip_case = True
                    break
            elif var_name in ds:
                data_month = ds[var_name].isel(Time=0).values

                if var_name in ['timeMonthly_avg_landIceFreshwaterFlux',
                                'timeMonthly_avg_landIceFreshwaterFluxTotal']:
                    sec_per_year = 365.25 * 24 * 3600
                    rho_fw = 1000.0
                    data_month = data_month * sec_per_year / rho_fw
                elif var_name == 'timeMonthly_avg_landIceFrictionVelocity':
                    m2cm = 100
                    data_month = data_month * m2cm

                monthly_data_arrays.append(data_month)
            else:
                print(f"--> Error: Variable '{var_name}' missing in dataset {os.path.basename(file_path)}.")
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
            dx=dx, cases_str=combined_cases_str, max_level_cell=max_level_cell,
            opt_region=opt_region, iis=iis, lon=lon, lat=lat,
            out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed, var_config=VAR_CONFIG
        )

    print(f"Finished processing and rendering all fields for: {dx}_{combined_cases_str} ({PLOT_VARIABLE}).\n")


# =========================================================================
# Helper function to dynamically retrieve configurations
# =========================================================================
def get_variable_config(variable_name, PLOT_ROSS=False):
    """Returns runtime mapping parameters based on global options string."""

    # 1. Assign the dictionary to a local variable instead of returning immediately
    if variable_name == 'Tbot':
        config = {
            'name': 'timeMonthly_avg_activeTracers_temperature',
            'vmin': -2.6, 'vmax': -1.6, 'contours': [], 'cmap': 'cmo.thermal',
            'cb_label': 'Sea Floor Temperature [°C]', 'title_prefix': 'Bottom Temperature',
            'file_prefix': 'Tbot', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Sbot':
        config = {
            'name': 'timeMonthly_avg_activeTracers_salinity',
            'vmin': 34.2, 'vmax': 35.0, 'contours': [], 'cmap': 'cmo.haline',
            'cb_label': 'Sea Floor Salinity [g/kg]', 'title_prefix': 'Bottom Salinity',
            'file_prefix': 'Sbot', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Tsurf':
        config = {
            'name': 'timeMonthly_avg_activeTracers_temperature',
            'vmin': -2.6, 'vmax': -1.6, 'contours': [], 'cmap': 'cmo.thermal',
            'cb_label': 'Surface Layer Temperature [°C]', 'title_prefix': 'Surface Temperature',
            'file_prefix': 'Tsurf', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Ssurf':
        config = {
            'name': 'timeMonthly_avg_activeTracers_salinity',
            'vmin': 34.2, 'vmax': 35.0, 'contours': [], 'cmap': 'cmo.haline',
            'cb_label': 'Surface Layer Salinity [g/kg]', 'title_prefix': 'Surface Salinity',
            'file_prefix': 'Ssurf', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'BotSpeed':
        config = {
            'name': 'BotSpeed',
            'vmin': 0.0, 'vmax': 0.3, 'contours': [], 'cmap': 'cmo.speed',
            'cb_label': 'Bottom Flow Speed [m/s]', 'title_prefix': 'Bottom Flow Speed',
            'file_prefix': 'BotSpeed', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'SurfSpeed':
        config = {
            'name': 'SurfSpeed',
            'vmin': 0.0, 'vmax': 0.3, 'contours': [], 'cmap': 'cmo.speed',
            'cb_label': 'Surface Flow Speed [m/s]', 'title_prefix': 'Surface Flow Speed',
            'file_prefix': 'SurfSpeed', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'DepthAvgSpeed':
        config = {
            'name': 'DepthAvgSpeed',
            'vmin': 0.0, 'vmax': 0.3, 'contours': [], 'cmap': 'cmo.speed',
            'cb_label': 'Depth Averaged Speed [m/s]', 'title_prefix': 'Depth Averaged Speed',
            'file_prefix': 'DepthAvgSpeed', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Tint':
        config = {
            'name': 'timeMonthly_avg_activeTracers_temperature',
            'vmin': -2.5, 'vmax': 0.0, 'contours': [], 'cmap': 'cmo.thermal',
            'cb_label': 'Depth Averaged Temperature [°C]', 'title_prefix': 'Depth Averaged Temperature',
            'file_prefix': 'Tint', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Sint':
        config = {
            'name': 'timeMonthly_avg_activeTracers_salinity',
            'vmin': 34.2, 'vmax': 35.0, 'contours': [], 'cmap': 'cmo.haline',
            'cb_label': 'Depth Averaged Salinity [g/kg]', 'title_prefix': 'Depth Averaged Salinity',
            'file_prefix': 'Sint', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'ColSpeed':
        config = {
            'name': 'timeMonthly_avg_columnIntegratedSpeed',
            'vmin': 0.0, 'vmax': 0.3, 'contours': [], 'cmap': 'cmo.speed',
            'cb_label': 'Column Integrated Speed [m/s]', 'title_prefix': 'Column Integrated Speed',
            'file_prefix': 'ColSpeed', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'MLD':
        config = {
            'name': 'timeMonthly_avg_dThreshMLD',
            'vmin': 0, 'vmax': 500, 'contours': [], 'cmap': 'cmo.deep',
            'cb_label': 'Mixed Layer Depth [m]', 'title_prefix': 'MLD',
            'file_prefix': 'MLD', 'opt_proj': 'wed', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'GMkappa':
        config = {
            'name': 'GMkappa',
            'vmin': 0.0, 'vmax': 1.0, 'contours': [], 'cmap': 'CMRmap_r',
            'cb_label': r'$\kappa_{GM}$ / max($\kappa_{GM}$)', 'title_prefix': r'Normalized $\kappa_{GM}$',
            'file_prefix': 'GMkappa', 'opt_proj': 'sps', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Melt':
        config = {
            'name': 'timeMonthly_avg_landIceFreshwaterFlux',
            'vmin': -3.0, 'vmax': 3.0, 'contours': [], 'cmap': 'RdBu_r',
            'cb_label': 'Melt rate interfacial [m/a]', 'title_prefix': 'Melt rate Interfacial',
            'file_prefix': 'Melt', 'opt_proj': 'fris', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'MeltTotal':
        config = {
            'name': 'timeMonthly_avg_landIceFreshwaterFluxTotal',
            'vmin': -3.0, 'vmax': 3.0, 'contours': [], 'cmap': 'RdBu_r',
            'cb_label': 'Melt rate total [m/a]', 'title_prefix': 'Melt rate total',
            'file_prefix': 'MeltTot', 'opt_proj': 'fris', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Ustar':
        config = {
            'name': 'timeMonthly_avg_landIceFrictionVelocity',
            'vmin': 0.4, 'vmax': 1.2, 'contours': [], 'cmap': 'cmo.speed',
            'cb_label': 'Ustar [cm/s]', 'title_prefix': 'Land Ice Friction Velocity (Ustar)',
            'file_prefix': 'Ustar', 'opt_proj': 'fris', 'vmax_scale_factor': 0.25
        }
    elif variable_name == 'Tstar':
        config = {
            'name': 'Tstar',
            'vmin': 0.0, 'vmax': 1.0, 'contours': [], 'cmap': 'hot_r',
            'cb_label': 'Boundary-Interface Temp Difference [°C]', 'title_prefix': 'Tstar Temperature Difference',
            'file_prefix': 'Tstar', 'opt_proj': 'fris', 'vmax_scale_factor': 0.1
        }
    else:
        raise ValueError(f"Unknown variable configuration requested: {variable_name}")

    # 2. Check the flag and overwrite the value if true
    if PLOT_ROSS:
        config['opt_proj'] = 'ross'

    return config


# =========================================================================
# 3. Main Execution Block
# =========================================================================
if __name__ == "__main__":

    fris_loc = '/pscratch/sd/v/vankova/fris_analysis/fris_plots/spatial_stats'
    opt_region = False
    iceshelves = ["Shelf"]

    # -----------------------------------------------------------------
    # SPECIFY TARGET SIMULATION TYPE AND YEARS
    # -----------------------------------------------------------------
    RUN_TYPE = 'Spin6'
    TARGET_YEARS = ['0002', '0003', '0004']

    # Array of target parameters to map out in parallel
    # PLOT_VARIABLES = ['Sbot', 'Sint', 'Tbot', 'Tint', 'ColSpeed', 'MLD', 'Tstar', 'Ustar', 'MeltTotal', 'Melt', 'Tsurf',
    #                  'Ssurf', 'BotSpeed', 'SurfSpeed', 'DepthAvgSpeed']

    # PLOT_VARIABLES = ['ColSpeed', 'MLD', 'Tstar', 'Ustar', 'MeltTotal', 'Melt']
    # PLOT_VARIABLES = ['Sbot', 'Sint', 'Tbot', 'Tint']
    PLOT_VARIABLES = ['Ssurf', 'Tsurf', 'SurfSpeed', 'BotSpeed', 'DepthAvgSpeed']

    PLOT_ROSS = False

    if PLOT_ROSS:
        fris_loc = f'{fris_loc}/ross'

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
    # Multi-dimensional bundling: Cross-product of simulations (resolutions) * chosen plot variables
    tasks = []
    for var in PLOT_VARIABLES:
        var_config = get_variable_config(var, PLOT_ROSS)
        for Fnum, cases in simulations.items():
            tasks.append(
                (Fnum, cases, var_config, RUN_TYPE, TARGET_YEARS, fris_loc, opt_region, iceshelves, var)
            )

    # Scale allocation max pool caps based on work list volume (usually 4 to 16 combinations)
    num_processes = min(len(tasks), 16)

    print(f"Spawning an isolated execution pool of {num_processes} parallel processes to process {len(tasks)} tasks...")

    with Pool(processes=num_processes) as pool:
        pool.map(process_single_resolution, tasks)

    print("All spatial resolution pools and multi-variable configurations completed successfully.")