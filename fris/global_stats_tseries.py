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


def parse_dt_to_seconds(dt_val):
    """Converts a config_dt value (string, byte, or number) to float seconds."""
    if dt_val is None:
        return None

    # Decode if it's a byte string
    if isinstance(dt_val, bytes):
        dt_val = dt_val.decode('utf-8')

    if isinstance(dt_val, str):
        dt_val = dt_val.strip().strip("'\"")
        # Handle hh:mm:ss format (e.g. "00:10:00")
        if ':' in dt_val:
            try:
                parts = dt_val.split(':')
                if len(parts) == 3:
                    h, m, s = map(float, parts)
                    return h * 3600.0 + m * 60.0 + s
                elif len(parts) == 2:
                    m, s = map(float, parts)
                    return m * 60.0 + s
            except ValueError:
                pass

        # If it's a plain numeric string (e.g. "600")
        try:
            return float(dt_val)
        except ValueError:
            print(f"WARNING: Could not parse string dt value: {dt_val}")
            return 0.0

    # If it's already a numeric type or a numpy/xarray scalar
    try:
        return float(dt_val)
    except (TypeError, ValueError):
        return 0.0


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
            if (idx + 1) % 10 == 0 or idx == 0 or (idx + 1) == len(file_list):
                print(f"[{idx + 1}/{len(file_list)}] Parsing: {file_path.split('/')[-1]}")

            with xr.open_dataset(file_path, decode_timedelta=True) as ds:
                file_metrics = {}

                # Drop spatial indexes; grab the global coordinate data cleanly
                for varname in vars_global:
                    if varname in ds:
                        file_metrics[varname] = ds[varname]
                    else:
                        alt_name = f"globalStats_{varname}"
                        if alt_name in ds:
                            file_metrics[varname] = ds[alt_name]

                # --- Extract config_dt and process it into a float ---
                raw_dt = None
                if 'config_dt' in ds.attrs:
                    raw_dt = ds.attrs['config_dt']
                elif 'config_dt' in ds:
                    raw_dt = ds['config_dt'].values

                # Parse raw string/byte format down to total numeric float seconds
                dt_seconds = parse_dt_to_seconds(raw_dt)

                if dt_seconds is not None:
                    ref_var = list(file_metrics.values())[0] if file_metrics else None
                    if ref_var is not None:
                        # Force creation of a clean float64 array matching the length of other fields
                        file_metrics['config_dt'] = xr.full_like(ref_var, dt_seconds, dtype=np.float64)
                else:
                    print(f"WARNING: 'config_dt' not found in {file_path.split('/')[-1]}")

                if file_metrics:
                    daily_ds = xr.Dataset(file_metrics)
                    daily_tseries_list.append(daily_ds)

        if daily_tseries_list:
            print("Concatenating daily global stats into master time-series dataset...")
            out_ds = xr.concat(daily_tseries_list, dim='Time')

            # Ensure 'Time' is the leading dimension
            if 'Time' in out_ds.dims:
                out_ds = out_ds.transpose('Time', ...)

            # -----------------------------------------------------------------
            # --- GAP HANDLING: Build a continuous, gap-aware Time Axis ---
            # -----------------------------------------------------------------
            try:
                # 1. Inspect the calendar type parsed by xarray (e.g., 'noleap', 'gregorian')
                if hasattr(out_ds['Time'].values[0], 'calendar'):
                    calendar = out_ds['Time'].values[0].calendar
                else:
                    calendar = 'noleap'  # Typical safe fallback for MPAS-O configurations

                # 2. Extract bounding limits directly from the file arrays
                start_time = out_ds['Time'].min().values.item()
                end_time = out_ds['Time'].max().values.item()

                # 3. Create a clean, gap-free target axis spanning the entire stretch.
                # Note: If your simulation steps are daily, swap 'MS' (Month Start) to 'D' (Day).
                full_time_axis = xr.cftime_range(
                    start=start_time,
                    end=end_time,
                    freq='MS',
                    calendar=calendar
                )

                # 4. Reindex alignment. Any missing steps automatically fill with NaN pointers.
                print(f"Reindexing timeline from {start_time} to {end_time} to expose data gaps...")
                out_ds = out_ds.reindex(Time=full_time_axis)

            except Exception as e:
                print(f"WARNING: Automatic timeline reindexing failed: {e}")
                print("Defaulting to non-reindexed array sequence.")
            # -----------------------------------------------------------------

            # Save dataset to NetCDF file
            output_filename = f'{out_prl}/global_stats_tseries_{dx}_{sec}{subsec}.nc'
            out_ds.to_netcdf(output_filename)
            print(f"Successfully consolidated global diagnostics to: {output_filename}")
        else:
            print(f"No valid variables extracted for configuration: {dx}_{sec}{subsec}")