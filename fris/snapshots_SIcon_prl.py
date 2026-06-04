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
def process_single_file(file_path, idx, total_files, dx, sec, subsec_str, opt_region, iis, lon, lat,
                        out_plot_dir, dsMesh_trimmed):
    """Processes a single mpassi history file and generates a sea ice concentration plot."""
    base_filename = os.path.basename(file_path)

    try:
        date_str = base_filename.split('.')[-2]
    except IndexError:
        date_str = f"idx_{idx:03d}"

    # Print with process ID to track progress in logs
    print(f"[{os.getpid()}] Processing {idx + 1}/{total_files}: {base_filename}")

    with xr.open_dataset(file_path) as ds:
        # Extract Sea Ice Concentration (iceAreaCell is the fraction between 0 and 1)
        if 'timeMonthly_avg_iceAreaCell' in ds:
            ice_conc = ds['timeMonthly_avg_iceAreaCell'].isel(Time=0).values
        elif 'iceAreaCell' in ds:
            ice_conc = ds['iceAreaCell'].isel(Time=0).values
        else:
            return f"Error: Sea ice concentration missing in {base_filename}"

    # Apply region restrictions if toggled
    if opt_region:
        ice_conc = np.where(iis == 0, np.nan, ice_conc)

    # Plot generation
    projection = ccrs.SouthPolarStereo()
    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True, subplot_kw={"projection": projection})

    descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)
    # Norm bounds mapped directly from 0 (open water) to 1 (solid pack ice)
    collection = mosaic.polypcolor(ax, descriptor, ice_conc, norm=colors.Normalize(vmin=0.0, vmax=1.0),
                                   cmap="cmo.ice")

    # Overlay standard 15% (0.15) Sea Ice Extent Outline Contour
    #valid_mask_ice = ~np.isnan(ice_conc)
    #if np.any(valid_mask_ice):
    #    ax.tricontour(
    #        lon[valid_mask_ice], lat[valid_mask_ice], ice_conc[valid_mask_ice],
    #        levels=[0.15], colors='black', linewidths=1.5, transform=ccrs.PlateCarree()
    #    )

    #ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

    fig.colorbar(collection, fraction=0.1, shrink=0.5, label="Sea Ice Area Fraction")
    ax.set_title(f"Sea Ice Concentration & 15% Extent Contour - {dx} {sec} {subsec_str} ({date_str})", fontsize=12,
                 pad=10)

    # Write output to sea ice snapshots directory
    output_png_path = f'{out_plot_dir}/IceConc_{dx}_{sec}{subsec_str}_{date_str}.png'
    plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
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

    simulations = {
        '8': [('Spin6', 'p1'), ('Spin1', 'p1')],
        '4': [('Spin6', 'p1'), ('Spin1', 'p1')],
        '2': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2')],
        '1': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
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

            # Resolution-separated snapshot location for Sea Ice
            out_plot_dir = f'{fris_loc}/snapshots_SIConc/{dx}_{sec}{subsec_str}'
            os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

            # SWAPPED PATTERN FROM mpaso TO mpassi
            file_pattern = f"{fpath}/{run_name}.mpassi.hist.am.timeSeriesStatsMonthly.*.nc"
            file_list = sorted(glob.glob(file_pattern))
            total_files = len(file_list)

            if not file_list:
                print(f"--> Warning: Missing sea ice files matching: {file_pattern}. Skipping.")
                continue

            print(f"\nSpinning up parallel pool with {NUM_WORKERS} workers for {total_files} sea ice files...")

            # Use partial to freeze non-iterable configuration arguments (removed max_level_cell)
            worker_task = partial(
                process_single_file,
                total_files=total_files, dx=dx, sec=sec, subsec_str=subsec_str,
                opt_region=opt_region, iis=iis, lon=lon, lat=lat,
                out_plot_dir=out_plot_dir, dsMesh_trimmed=dsMesh_trimmed
            )

            # Map the file list across the process pool
            iterable_args = [(file_path, idx) for idx, file_path in enumerate(file_list)]

            with mp.Pool(processes=NUM_WORKERS) as pool:
                results = pool.starmap(worker_task, iterable_args)

            # Print any individual worker errors
            for res in results:
                if res is not None:
                    print(res)

            print(f"Completed rendering all sea ice frames for configuration: {dx}_{sec}{subsec_str}.")