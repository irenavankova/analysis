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

# Define global-averaged diagnostic variables requested
vars_global = [
    'CFLNumberGlobal',
    'kineticEnergyCellMax'
]


# -------------------------------------------------------------------------
# Preprocess function for high-performance lazy loading
# -------------------------------------------------------------------------
def extract_vars(ds):
    """
    Slices out only the requested global metrics immediately upon opening
    each file. This avoids choking memory when opening thousands of files.
    """
    # Detect matching variables or alternatives with 'globalStats_' prefixes
    valid_vars = [v for v in vars_global if v in ds]
    alt_vars = [f"globalStats_{v}" for v in vars_global if f"globalStats_{v}" in ds]

    # Extract dataset subset containing ONLY the target variables
    ds_subset = ds[valid_vars + alt_vars]

    # Normalize naming convention if prefixes were applied by the AM module
    rename_dict = {f"globalStats_{v}": v for v in vars_global if f"globalStats_{v}" in ds_subset}
    if rename_dict:
        ds_subset = ds_subset.rename(rename_dict)

    return ds_subset


# -------------------------------------------------------------------------
# Main Execution Loop
# -------------------------------------------------------------------------
for Fnum, cases in simulations.items():
    dx = f'F{Fnum}'

    for sec, subsec in cases:
        print("\n" + "=" * 60)
        print(f"STARTING CONFIGURATION: Resolution={dx} | Section={sec} | Sub-section={subsec}")
        print("=" * 60)

        # Build execution paths identically to your original setup
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
            print(f"WARNING: No global stats files found for pattern: {file_pattern}. Skipping.")
            continue

        print(f"Found {len(file_list)} files. Opening dataset and processing in parallel...")

        # --- Step 3: Fast Parallel Metadata Aggregation ---
        # --- Step 3: Fast Parallel Metadata Aggregation ---
        try:
            # open_mfdataset reads NetCDF headers concurrently using local Dask worker threads
            # Added decode_timedelta=True to silence the Python 3.14 warnings
            out_ds = xr.open_mfdataset(
                file_list,
                concat_dim='Time',
                combine='nested',
                preprocess=extract_vars,
                parallel=True,
                data_vars='minimal',
                coords='minimal',
                compat='override',
                decode_timedelta=True  # <--- Silences FutureWarnings
            )

            # Ensure proper dimension ordering (Time as leading axis)
            if 'Time' in out_ds.dims:
                out_ds = out_ds.transpose('Time', ...)

            # Define the destination path
            output_filename = f'{out_prl}/global_stats_tseries_{dx}_{sec}{subsec}.nc'

            # CRITICAL FIX: Load data into memory first to bypass Dask multithreading
            # conflicts inside the HDF5/NetCDF4 engine during the write phase.
            print("Loading aggregated timeseries into memory...")
            out_ds.load()

            # Compute and save consolidated time series data cleanly
            print(f"Writing aggregated timeseries out to netCDF...")
            out_ds.to_netcdf(output_filename)
            print(f"SUCCESS: Consolidated global diagnostics saved to: {output_filename}")

            # Explicitly close dataset to release system file handles
            out_ds.close()

        except Exception as e:
            print(f"ERROR: Failed processing configuration {dx}_{sec}{subsec}: {e}")