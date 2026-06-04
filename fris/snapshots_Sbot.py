#!/usr/bin/env python3

import glob
import os
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean
import mosaic
import matplotlib.path as mpath
from matplotlib import colors

# =========================================================================
# 1. Configuration & Global Path Frameworks
# =========================================================================
fris_loc = '/Users/ivankova/Desktop/Fris_hr'
opt_region = False
iceshelves = ["Shelf"]

# Dictionary structuring resolutions and simulation phases (from obs_tseries.py)
simulations = {
    '8': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '4': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '2': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2')],
    '1': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
}

# =========================================================================
# 2. Main Simulation & Resolution Processing Loop
# =========================================================================
for Fnum, cases in simulations.items():
    dx = f'F{Fnum}'

    # ---------------------------------------------------------------------
    # DYNAMIC MESH GENERATION (Mirrors obs_tseries.py logic)
    # ---------------------------------------------------------------------
    # Construct the exact path to the Spin6 restart file to read the grid mesh
    run_name_mask = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
    fpath_mask = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name_mask}/run'
    mesh_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

    if not os.path.exists(mesh_file):
        print(f"--> Warning: Mesh file missing for {dx}: {mesh_file}. Skipping resolution.")
        continue

    print(f"\nParsing Dynamic Mesh Properties for Resolution {dx} via Spin6 Restart Grid...")
    with xr.open_dataset(mesh_file) as dsMesh:
        # Load coordinate geometry upfront to reduce disk IO bottlenecks
        dsMesh_trimmed = dsMesh[['latCell', 'lonCell', 'landIceFloatingMask', 'areaCell',
                                 'maxLevelCell', 'nEdges', 'nVertices', 'lonEdge', 'latEdge',
                                 'lonVertex', 'latVertex', 'cellsOnEdge', 'cellsOnVertex',
                                 'verticesOnEdge', 'verticesOnCell', 'edgesOnVertex']].load()

    # Geometry Extractions
    lat = dsMesh_trimmed['latCell'].values * 180 / np.pi
    lon = dsMesh_trimmed['lonCell'].values * 180 / np.pi
    max_level_cell = dsMesh_trimmed['maxLevelCell'].values - 1  # Convert 1-indexed to 0-indexed

    # Apply region restrictions if toggled
    if opt_region:
        import gmask_reg

        iam = gmask_reg.get_mask(iceshelves, mesh_file, opt_noGL=0, opt_wct=1)
        iis = iam[0, :]

    # ---------------------------------------------------------------------
    # Configuration Phases (Spin6, Spin1 sub-cases)
    # ---------------------------------------------------------------------
    for sec, subsec in cases:
        subsec_str = subsec if sec == 'Spin1' else ''
        print("\n" + "=" * 60)
        print(f"PROCESSING RUN: Res={dx} | Section={sec} | Sub-section={subsec}")
        print("=" * 60)

        # Dynamic root path allocation (from obs_tseries.py)
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

        # Set up output directories specific to this resolution/simulation subset
        out_plot_dir = f'{fris_loc}/Fris_plots/snapshots_xy/{dx}_{sec}{subsec_str}'
        os.makedirs(out_plot_dir, mode=0o755, exist_ok=True)

        # Glob history archives
        file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"
        file_list = sorted(glob.glob(file_pattern))

        if not file_list:
            print(f"--> Warning: Missing files matching: {file_pattern}. Skipping phase.")
            continue

        print(f"Found {len(file_list)} target history files. Beginning mapping loop...")

        # =========================================================================
        # 3. History File Parsing & Plot Compilation
        # =========================================================================
        for idx, file_path in enumerate(file_list):
            base_filename = os.path.basename(file_path)

            # Extract timestamp string (e.g., "0002-01-01") for output naming
            try:
                date_str = base_filename.split('.')[-2]
            except IndexError:
                date_str = f"idx_{idx:03d}"

            print(f"   [{idx + 1}/{len(file_list)}] Rendering Bottom Salinity for: {base_filename}")

            with xr.open_dataset(file_path) as ds:
                # Isolate target key strings safely
                if 'timeMonthly_avg_activeTracers_salinity' in ds:
                    salinity_3d = ds['timeMonthly_avg_activeTracers_salinity'].isel(Time=0).values
                elif 'salinity' in ds:
                    salinity_3d = ds['salinity'].isel(Time=0).values
                else:
                    print(f"   --> Error: Salinity key configuration mismatch in {base_filename}. Skipping file.")
                    continue

            # Advanced numpy indexing to slice the bottom layer for every mesh cell
            num_cells = salinity_3d.shape[0]
            bottom_salinity = salinity_3d[np.arange(num_cells), max_level_cell]

            # Mask uninitialized / invalid cells out
            bottom_salinity = np.where(max_level_cell >= 0, bottom_salinity, np.nan)

            if opt_region:
                bottom_salinity = np.where(iis == 0, np.nan, bottom_salinity)

            # Map Limits and Colormap Specifications
            v_min, v_max = 34.0, 35.0
            norm = colors.Normalize(vmin=v_min, vmax=v_max)
            varcmap = "cmo.haline"
            plotlabel = "Sea Floor Salinity [psu]"

            # Map Generation Initializations
            projection = ccrs.PlateCarree()
            fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True, subplot_kw={"projection": projection})

            # Native MPAS unstructured mesh interpolation via Mosaic
            descriptor = mosaic.Descriptor(dsMesh_trimmed, projection=projection)
            collection = mosaic.polypcolor(ax, descriptor, bottom_salinity, norm=norm, cmap=varcmap)

            # Overlay Specific Isohaline Contours (e.g., 34.9 psu critical threshold)
            valid_mask = ~np.isnan(bottom_salinity)
            if np.any(valid_mask):
                ax.tricontour(
                    lon[valid_mask], lat[valid_mask], bottom_salinity[valid_mask],
                    levels=[34.9], colors='black', linewidths=1.5, transform=ccrs.PlateCarree()
                )

            # Establish Spatial Domain Extents & Aesthetics
            ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())
            ax.set_aspect('auto')
            ax.gridlines(draw_labels=True)
            ax.set_facecolor('lightgray')

            # Incorporate Colorbar Metrics & Plot Header Title
            fig.colorbar(collection, fraction=0.1, shrink=0.5, label=plotlabel)
            ax.set_title(f"Bottom Salinity - {dx} {sec} {subsec_str} ({date_str})", fontsize=12, pad=10)

            # Save clean resolution-separated images
            output_png_path = f'{out_plot_dir}/Sbot_{dx}_{sec}{subsec_str}_{date_str}.png'
            plt.savefig(output_png_path, bbox_inches='tight', dpi=200)
            plt.close(fig)  # Clear system memory references

        print(f"Completed rendering all frames for configuration: {dx}_{sec}{subsec_str}.")