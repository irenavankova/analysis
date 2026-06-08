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
        # 1. Extract the current calendar type (e.g., 'noleap', 'gregorian')
        # If Xarray decoded it as an object array of cftime objects:
        if hasattr(out_ds['Time'].values[0], 'calendar'):
            calendar = out_ds['Time'].values[0].calendar
        else:
            calendar = 'noleap'  # Safe default fallback for MPAS

        # 2. Find the strict start and end dates present in your files
        start_time = out_ds['Time'].min().values.item()
        end_time = out_ds['Time'].max().values.item()

        # 3. Create a continuous, gap-free time axis at daily ('D') or monthly ('MS') frequency
        # Given your file example ends with '0004-06-01', '0004-07-01', it looks like monthly steps:
        # If they are daily, change 'MS' (Month Start) to 'D' (Day)
        full_time_axis = xr.cftime_range(
            start=start_time,
            end=end_time,
            freq='MS',  # Use 'D' if your data is daily instead of monthly
            calendar=calendar
        )

        # 4. Reindex the dataset. Gaps will automatically fill with NaN
        print(f"Reindexing timeline from {start_time} to {end_time} to expose gaps...")
        out_ds = out_ds.reindex(Time=full_time_axis)

    except Exception as e:
        print(f"WARNING: Could not automatically reindex time axis due to: {e}")
        print("Proceeding with concatenated array (gaps will be skipped rather than NaN-filled).")
    # -----------------------------------------------------------------

    # Save dataset to NetCDF file
    output_filename = f'{out_prl}/global_stats_tseries_{dx}_{sec}{subsec}.nc'
    out_ds.to_netcdf(output_filename)
    print(f"Successfully consolidated global diagnostics to: {output_filename}")
else:
    print(f"No valid variables extracted for configuration: {dx}_{sec}{subsec}")