#!/usr/bin/env python3

import glob
import xarray as xr
import numpy as np
import gmask_reg

# ==============================================================================
# CONFIGURATION OPTIONS
# ==============================================================================
# Define the structured matrix of runs to loop through
# Format: { Fnum: [(sec, subsec), ...] }
simulations = {
    '8': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '4': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '2': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2')],
    '1': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
}

opt_noGL = 1
exclude_high_vol_cells = False  # Set to True to exclude cells with volume > 100, False to keep all cells

units_scale_factor = 60 * 60 * 24 * 365
regs = ["FRISshelf", "RonneDshelf", "FilchnerDshelf", "BerknerBank"]

# ==============================================================================
# MASTER RUN LOOP
# ==============================================================================
for Fnum, cases in simulations.items():
    dx = f'F{Fnum}'

    for sec, subsec in cases:
        print("\n" + "=" * 60)
        print(f"STARTING CONFIGURATION: Resolution={dx} | Section={sec} | Sub-section={subsec}")
        print(f"Filtering Out High-Volume Cells (>100): {exclude_high_vol_cells}")
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

        fpath_mask = fpath
        run_name_mask = run_name

        # Define file lists
        file_pattern = f"{fpath}/{run_name}.mpassi.hist.am.timeSeriesStatsMonthly.*.nc"
        file_list = sorted(glob.glob(file_pattern))

        if not file_list:
            print(f"Warning: No history files found matching pattern: {file_pattern}. Skipping simulation case...")
            continue

        # Step 1: Load the mask and grid metrics from the ocean restart file
        mask_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

        try:
            mask_ds = xr.open_dataset(mask_file)
        except FileNotFoundError:
            print(f"Warning: Restart file {mask_file} not found. Skipping simulation case...")
            continue

        areaCell = mask_ds['areaCell']  # Horizontal grid cell area (nCells,)

        print(f"Getting masks...")
        iam = gmask_reg.get_mask(regs, mask_file, opt_noGL=opt_noGL)

        print(f"Reshaping masks...")
        iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': regs})

        # --- Dynamic Pre-screening Choice ---
        if exclude_high_vol_cells:
            print("Pre-screening files to find cells with sea-ice volume > 100 at any point...")
            cells_to_exclude = np.zeros(areaCell.shape, dtype=bool)

            for file_path in file_list:
                with xr.open_dataset(file_path, decode_timedelta=True) as ds_check:
                    vol_check = ds_check['timeMonthly_avg_iceVolumeCell']
                    high_vol_mask = (vol_check > 100.0).any(dim='Time').values
                    cells_to_exclude = cells_to_exclude | high_vol_mask

            num_excluded = np.sum(cells_to_exclude)
            print(f"Excluding {num_excluded} cells out of {len(cells_to_exclude)} due to the volume > 100 threshold.")
            valid_cells_mask = xr.DataArray(~cells_to_exclude, dims=['nCells'])
        else:
            print("Skipping pre-screening loop. All spatial cells will be included in calculations.")
            valid_cells_mask = xr.DataArray(np.ones(areaCell.shape, dtype=bool), dims=['nCells'])

        # List to accumulate the lightweight time series results from each file
        monthly_tseries_list = []

        print(f"Found {len(file_list)} files. Starting sequential processing...")

        # Step 2: Loop over files one by one to keep memory profile flat
        for idx, file_path in enumerate(file_list):
            print(f"[{idx + 1}/{len(file_list)}] Processing: {file_path.split('/')[-1]}")

            with xr.open_dataset(file_path, decode_timedelta=True) as ds:

                # 1. Base State Variables
                ice_concentration = ds['timeMonthly_avg_iceAreaCell']
                ice_thickness = ds['timeMonthly_avg_iceVolumeCell']

                # 2. Extract and Calculate production contributing terms (Scaled to m/yr)
                congelation = ds['timeMonthly_avg_congelation'] * units_scale_factor
                frazil = ds['timeMonthly_avg_frazilFormation'] * units_scale_factor
                snowice = ds['timeMonthly_avg_snowiceFormation'] * units_scale_factor

                # Combined Production
                ice_production = congelation + frazil + snowice

                # 3. Extract and Calculate Total Sea Ice Melting
                basal = ds['timeMonthly_avg_basalIceMelt'] * units_scale_factor
                surface = ds['timeMonthly_avg_surfaceIceMelt'] * units_scale_factor
                lateral = ds['timeMonthly_avg_lateralIceMelt'] * units_scale_factor
                ice_melting = basal + surface + lateral

                file_metrics = {}

                # Dictionary of fields where we want mean, max, and min extraction
                vars_to_process = {
                    'ice_concentration': ice_concentration,
                    'ice_thickness': ice_thickness,
                    'ice_production': ice_production,
                    'ice_melting': ice_melting,
                    'congelation': congelation,
                    'frazil': frazil,
                    'snowice': snowice
                }

                # Broadcast grid area safely across regions and apply masks once per file loop
                area_broadcasted = areaCell.expand_dims(region=iam_da.region)
                area_masked = area_broadcasted.where(iam_da).where(valid_cells_mask)

                # ---- Calculate Mean, Max, Min dynamically for all variables ----
                for name, var in vars_to_process.items():
                    var_broadcasted = var.expand_dims(region=iam_da.region)
                    var_masked = var_broadcasted.where(iam_da).where(valid_cells_mask)

                    # A. Mean (Area-Weighted)
                    active_area = area_masked.where(var_masked.notnull())
                    weighted_sum = (var_masked * active_area).sum(dim='nCells', skipna=True)
                    total_area = active_area.sum(dim='nCells', skipna=True)
                    file_metrics[f'mean_{name}'] = (weighted_sum / total_area).transpose('Time', 'region')

                    # B. Maximum
                    file_metrics[f'max_{name}'] = var_masked.max(dim='nCells', skipna=True).transpose('Time', 'region')

                    # C. Minimum
                    file_metrics[f'min_{name}'] = var_masked.min(dim='nCells', skipna=True).transpose('Time', 'region')

                # ---- Calculate Integrated Ice Volume ----
                thick_broadcasted = ice_thickness.expand_dims(region=iam_da.region)
                thick_masked = thick_broadcasted.where(iam_da).where(valid_cells_mask)
                file_metrics['integrated_ice_volume'] = (thick_masked * area_masked).sum(dim='nCells',
                                                                                         skipna=True).transpose('Time',
                                                                                                                'region')

                # Convert the dictionary of metrics into a single dataset and append to accumulator
                monthly_ds = xr.Dataset(file_metrics)
                monthly_tseries_list.append(monthly_ds)

        print("Finished processing individual files. Concatenating time series...")

        # Step 3: Concatenate all the processed monthly intervals along the 'Time' dimension
        out_ds = xr.concat(monthly_tseries_list, dim='Time')

        # Ensure correct coordinate order
        out_ds = out_ds.transpose('Time', 'region')

        # ---- ADD METADATA ATTRIBUTES ----
        # Assign config options as NetCDF global attributes
        out_ds.attrs['opt_noGL'] = int(opt_noGL)
        out_ds.attrs['exclude_high_vol_cells'] = str(exclude_high_vol_cells)
        out_ds.attrs['high_volume_threshold'] = 100.0 if exclude_high_vol_cells else "N/A"

        # Save cleaner dataset to NetCDF
        output_filename = f'bulk_seaice_tseries_{dx}_{sec}{subsec_str}.nc'
        out_ds.to_netcdf(output_filename)
        print(f"Successfully saved consolidated time series to: {output_filename}")