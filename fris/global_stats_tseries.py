#!/usr/bin/env python3
import glob
import os
import xarray as xr
import numpy as np

# Output directory path
out_prl = '/global/homes/v/vankova/data_analysis/my_scripts/nc_files/global_stats'

# -------------------------------------------------------------------------
# Ensure output directory exists with explicit permissions
# -------------------------------------------------------------------------
desired_mode = 0o755

if not os.path.exists(out_prl):
    os.makedirs(out_prl, mode=desired_mode, exist_ok=True)
    print(f"Created directory: {out_prl} with permissions {oct(desired_mode)}")
else:
    os.chmod(out_prl, desired_mode)
    print(f"Directory already exists. Updated permissions to {oct(desired_mode)}")

# Keep the same structural tracking for resolution variants and branches
simulations = {
    '8': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '4': [('Spin6', 'p1'), ('Spin1', 'p1')],
    '2': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2')],
    '1': [('Spin6', 'p1'), ('Spin1', 'p1'), ('Spin1', 'p2'), ('Spin1', 'p3')]
}

for Fnum, cases in simulations.items():
    dx = f'F{Fnum}'

    for sec, subsec in cases:
        print("\n" + "=" * 60)
        print(f"STARTING CONFIGURATION: Resolution={dx} | Section={sec} | Sub-section={subsec}")
        print("=" * 60)

        # Build paths identically to your setup
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

        # --- Step 2: Target the globalStats historical daily files ---
        print("Scanning daily global stats files...")
        file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.globalStats.*.nc"
        file_list = sorted(glob.glob(file_pattern))

        if not file_list:
            print(f"WARNING: No global stats files found for: {file_pattern}. Skipping.")
            continue

        # Global-averaged diagnostic variables requested
        vars_global = [
            'CFLNumberGlobal',
            'kineticEnergyCellMax'
        ]

        daily_tseries_list = []
        print(f"Found {len(file_list)} files. Extracting global metrics...")

        # --- Step 3: Parse metrics sequentially ---
        for idx, file_path in enumerate(file_list):
            # Print update every 10 files to keep log cleaner since daily files are numerous
            if (idx + 1) % 10 == 0 or idx == 0 or (idx + 1) == len(file_list):
                print(f"[{idx + 1}/{len(file_list)}] Parsing: {file_path.split('/')[-1]}")

            with xr.open_dataset(file_path, decode_timedelta=True) as ds:
                file_metrics = {}

                # Drop spatial indexes; grab the global coordinate data cleanly
                for varname in vars_global:
                    if varname in ds:
                        file_metrics[varname] = ds[varname]
                    else:
                        # Fallback case if variable naming contains prefixes in the AM module
                        alt_name = f"globalStats_{varname}"
                        if alt_name in ds:
                            file_metrics[varname] = ds[alt_name]

                if file_metrics:
                    daily_ds = xr.Dataset(file_metrics)
                    daily_tseries_list.append(daily_ds)

        if daily_tseries_list:
            print("Concatenating daily global stats into master time-series dataset...")
            out_ds = xr.concat(daily_tseries_list, dim='Time')

            # Ensure 'Time' is the leading dimension
            if 'Time' in out_ds.dims:
                out_ds = out_ds.transpose('Time', ...)

            # Save dataset to NetCDF file
            output_filename = f'{out_prl}/global_stats_tseries_{dx}_{sec}{subsec}.nc'
            out_ds.to_netcdf(output_filename)
            print(f"Successfully consolidated global diagnostics to: {output_filename}")
        else:
            print(f"No valid variables extracted for configuration: {dx}_{sec}{subsec}")