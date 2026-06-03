#!/usr/bin/env python3
import glob
import os  # <-- Add this import
import xarray as xr
import numpy as np

# Import the coordinates directly from your existing script
from fris_coordinates import sites_config, find_nearest_mpas_cells, load_transect_as_sites

opt_sites = 'filchdepr'
loc_transects = '/global/homes/v/vankova/data_analysis/files'
out_prl = f'/global/homes/v/vankova/data_analysis/my_scripts/nc_files/pts_{opt_sites}'

# -------------------------------------------------------------------------
# Ensure out_prl directory exists with explicit permissions
# -------------------------------------------------------------------------
# Define the numeric octal permissions you want (e.g., 0o755 or 0o770)
desired_mode = 0o755

if not os.path.exists(out_prl):
    # Create directory with permissions. exist_ok=True prevents race condition errors
    os.makedirs(out_prl, mode=desired_mode, exist_ok=True)
    print(f"Created directory: {out_prl} with permissions {oct(desired_mode)}")
else:
    # If it already exists, explicitly enforce/update the permissions
    os.chmod(out_prl, desired_mode)
    print(f"Directory already exists. Updated permissions to {oct(desired_mode)}")

if opt_sites == 'obs':
    sites_extract = sites_config
elif opt_sites == 'shelfbreak':
    transect_nc = f'{loc_transects}/{opt_sites}.nc'
    sites_extract = load_transect_as_sites(transect_nc, name_prefix="SB", site_type = "cherry")
elif opt_sites == 'filchdepr':
    transect_nc = f'{loc_transects}/{opt_sites}.nc'
    sites_extract = load_transect_as_sites(transect_nc, name_prefix="FD", site_type = "cherry")
elif opt_sites == 'ronnedepr':
    transect_nc = f'{loc_transects}/{opt_sites}.nc'
    sites_extract = load_transect_as_sites(transect_nc, name_prefix="RD", site_type="cherry")
elif opt_sites == 'ronnecenter':
    transect_nc = f'{loc_transects}/{opt_sites}.nc'
    sites_extract = load_transect_as_sites(transect_nc, name_prefix="RC", site_type="cherry")
elif opt_sites == 'berknerwest':
    transect_nc = f'{loc_transects}/{opt_sites}.nc'
    sites_extract = load_transect_as_sites(transect_nc, name_prefix="BW", site_type="cherry")


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
    mask_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

    with xr.open_dataset(mask_file) as dsMesh:
        site_da = find_nearest_mpas_cells(sites_extract, dsMesh)

    for sec, subsec in cases:
        print("\n" + "=" * 60)
        print(f"STARTING CONFIGURATION: Resolution={dx} | Section={sec} | Sub-section={subsec}")
        print("=" * 60)

        # Default initialization for Spin6
        if sec == 'Spin6':
            subsec_str = ''
            run_name = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
            fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'

        elif sec == 'Spin1':
            subsec_str = subsec
            if Fnum == '8':
                run_name = "20231114.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC08to60E3r1.spinup.chicoma-cpu"
            elif Fnum == '4':
                run_name = "20231108.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC04to60E3r1.spinup.chicoma-cpu"
            elif Fnum == '2':
                if subsec == 'p1':
                    run_name = "20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.chicoma-cpu"
                else:  # p2
                    run_name = "20231208.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.anvil"
            elif Fnum == '1':
                if subsec == 'p1':
                    run_name = "20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.chicoma-cpu"
                elif subsec == 'p2':
                    run_name = "20231209.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.anvil"
                else:  # p3
                    run_name = "20240201.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinupY5.chicoma-cpu"

            fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY1/{run_name}/run'

        # --- Step 2: Set up historical variables and file targets ---
        print("Scanning history files...")
        file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"
        file_list = sorted(glob.glob(file_pattern))

        if not file_list:
            raise FileNotFoundError(f"No history files found matching pattern: {file_pattern}")

        # Separate target fields by structural dimensionality
        vars_3d = [
            'timeMonthly_avg_layerThickness',
            'timeMonthly_avg_velocityMeridional',
            'timeMonthly_avg_velocityZonal',
            'timeMonthly_avg_activeTracers_temperature',
            'timeMonthly_avg_activeTracers_salinity',
            'timeMonthly_avg_potentialDensity',
            'timeMonthly_avg_BruntVaisalaFreqTop',
            'timeMonthly_avg_activeTracersTend_temperatureTend',
            'timeMonthly_avg_activeTracersTend_salinityTend',
            'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_temperatureHorizontalAdvectionTendency',
            'timeMonthly_avg_activeTracerHorizontalAdvectionTendency_salinityHorizontalAdvectionTendency',
            'timeMonthly_avg_activeTracerVerticalAdvectionTendency_temperatureVerticalAdvectionTendency',
            'timeMonthly_avg_activeTracerVerticalAdvectionTendency_salinityVerticalAdvectionTendency',
            'timeMonthly_avg_activeTracerVertMixTendency_temperatureVertMixTendency',
            'timeMonthly_avg_activeTracerVertMixTendency_salinityVertMixTendency',
            'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency',
            'timeMonthly_avg_activeTracerHorMixTendency_salinityHorMixTendency',
            'timeMonthly_avg_activeTracerSurfaceFluxTendency_temperatureSurfaceFluxTendency',
            'timeMonthly_avg_activeTracerSurfaceFluxTendency_salinitySurfaceFluxTendency',
            'timeMonthly_avg_temperatureShortWaveTendency',
            'timeMonthly_avg_activeTracerNonLocalTendency_temperatureNonLocalTendency',
            'timeMonthly_avg_activeTracerNonLocalTendency_salinityNonLocalTendency',
            'timeMonthly_avg_activeTracerVerticalAdvectionTopFlux_temperatureVerticalAdvectionTopFlux',
            'timeMonthly_avg_activeTracerVerticalAdvectionTopFlux_salinityVerticalAdvectionTopFlux'
        ]

        # NEW: Added variables utilizing the interface vertical dimension (nVertLevelsP1)
        vars_3dp1 = [
            'timeMonthly_avg_vertVelocityTop',
            'timeMonthly_avg_vertMLEBolusVelocityTop',
            'timeMonthly_avg_vertGMBolusVelocityTop',
            'timeMonthly_avg_vertDiffTopOfCell',
            'timeMonthly_avg_vertViscTopOfCell'
        ]

        vars_2d = [
            'timeMonthly_avg_landIceFreshwaterFlux',
            'timeMonthly_avg_landIceFreshwaterFluxTotal',
            'timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature',
            'timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature',
            'timeMonthly_avg_landIceFrictionVelocity',
            'timeMonthly_avg_ssh',
            'timeMonthly_avg_dThreshMLD',
            'timeMonthly_avg_tThreshMLD',
            'timeMonthly_avg_boundaryLayerDepth'
        ]

        monthly_tseries_list = []
        print(f"Found {len(file_list)} files. Beginning incremental multi-point extraction...")

        # --- Step 3: Stream and isolate localized profiles sequentially ---
        for idx, file_path in enumerate(file_list):
            print(f"[{idx + 1}/{len(file_list)}] Parsing: {file_path.split('/')[-1]}")

            with xr.open_dataset(file_path, decode_timedelta=True) as ds:
                file_metrics = {}

                # Pull 2D fields cleanly at specified location indices
                for varname in vars_2d:
                    if varname in ds:
                        file_metrics[varname] = ds[varname].isel(nCells=site_da)

                # Pull 3D mid-layer fields safely
                for varname in vars_3d:
                    if varname in ds:
                        file_metrics[varname] = ds[varname].isel(nCells=site_da)

                # NEW: Pull 3D interface-layer fields safely (nVertLevelsP1)
                for varname in vars_3dp1:
                    if varname in ds:
                        file_metrics[varname] = ds[varname].isel(nCells=site_da)

                # Bundle file metrics together, preserving native time frames
                monthly_ds = xr.Dataset(file_metrics)
                monthly_tseries_list.append(monthly_ds)

        print("Concatenating monthly site profiles into master dataset...")
        out_ds = xr.concat(monthly_tseries_list, dim='Time')

        # Fix: Order dimensions logically, accounting for both nVertLevels and nVertLevelsP1
        # Variables will be ordered (Time, site, nVertLevels) OR (Time, site, nVertLevelsP1) based on their identity.
        final_dim_order = [d for d in ['Time', 'site', 'nVertLevels', 'nVertLevelsP1'] if d in out_ds.dims]
        out_ds = out_ds.transpose(*final_dim_order)

        # Save cleaner dataset to NetCDF
        output_filename = f'{out_prl}/pts_{opt_sites}_tseries_{dx}_{sec}{subsec}.nc'
        out_ds.to_netcdf(output_filename)
        print(f"Successfully consolidated point profiles to: {output_filename}")