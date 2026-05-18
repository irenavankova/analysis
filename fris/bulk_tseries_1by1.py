#!/usr/bin/env python3

import glob
import xarray as xr
import numpy as np

import gmask_reg

Fnum = '8'
dx = f'F{Fnum}'
sec = 'Spin1'
subsec = 'p1'



run_name = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'
fpath_mask = fpath
run_name_mask = run_name
if sec == 'Spin6':
    subsec = ''


if sec == 'Spin1':
    if Fnum == '8':
        run_name = f"20231114.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC08to60E3r1.spinup.chicoma-cpu"
        subsec = ''
    elif Fnum == '4':
        run_name = f"20231108.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC04to60E3r1.spinup.chicoma-cpu"
        subsec = ''
    elif Fnum == '2':
        if subsec == 'p1':
            run_name = f"20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.chicoma-cpu"
        else:
            run_name = f"20231208.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC02to60E3r1.spinup.anvil"
    elif Fnum == '1':
        if subsec == 'p1':
            run_name = f"20231118.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.chicoma-cpu"
        elif subsec == 'p2':
            run_name = f"20231209.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinup.anvil"
        else:
            run_name = f"20240201.GMPAS-JRA1p5-DIB-PISMF-TMIX.TL319_FRISwISC01to60E3r1.spinupY5.chicoma-cpu"

    fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY1/{run_name}/run'


# Step 1: Load the mask and grid metrics from the restart NetCDF file
mask_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

mask_ds = xr.open_dataset(mask_file)
areaCell = mask_ds['areaCell']  # Horizontal grid cell area (nCells,)

#regs = ["RonneDcavity", "FilchnerDcavity", "BerknerSouth"]
regs = ["FRIS", "FRISshelf", "RonneDshelf", "FilchnerDshelf", "BerknerBank", "RonneDcavity", "FilchnerDcavity", "BerknerSouth"]

print(f"Getting masks...")

iam = gmask_reg.get_mask(regs, mask_file)

print(f"Reshaping masks...")

# Reshape the region mask into an xarray DataArray for seamless broadcasting: (region, nCells)
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': regs})

print(f"Getting a file list...")

# Get a sorted list of all individual history files
#file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.0004.*.nc"
file_pattern = f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc"

file_list = sorted(glob.glob(file_pattern))

if not file_list:
    raise FileNotFoundError(f"No history files found matching pattern: {file_pattern}")

# List to accumulate the lightweight time series results from each file
monthly_tseries_list = []

print(f"Found {len(file_list)} files. Starting sequential processing...")

# Step 2: Loop over files one by one to keep memory profile flat
for idx, file_path in enumerate(file_list):
    print(f"[{idx + 1}/{len(file_list)}] Processing: {file_path.split('/')[-1]}")

    # Use 'with' context manager to guarantee the file handles close and memory is freed immediately
    with xr.open_dataset(file_path, decode_timedelta=True) as ds:

        # Pull variables explicitly for this single timestep
        T_3d = ds['timeMonthly_avg_activeTracers_temperature']
        S_3d = ds['timeMonthly_avg_activeTracers_salinity']
        Rho_3d = ds['timeMonthly_avg_potentialDensity']
        h_3d = ds['timeMonthly_avg_layerThickness']

        # Calculate 3D Cell Volume for this specific timestep
        cell_volume_3d = areaCell * h_3d

        # Local dictionary to collect variables from this file
        file_metrics = {}

        vars_3d = {
            'temp': T_3d,
            'salt': S_3d,
            'rho': Rho_3d
        }

        for name, var in vars_3d.items():
            # 1. Broadcast variable across horizontal region mask dimensions
            print(f"Broadcasting regions...")
            var_broadcasted = var.expand_dims(region=iam_da.region)

            # 2. Mask out data points falling outside boundaries or that are uninitialized
            var_masked = var_broadcasted.where(iam_da).where(var_broadcasted.notnull())

            print(f"Calculating metrics...")
            # 3. Calculate spatial Minimum and Maximum over both spatial dimensions (nCells, nVertLevels)
            file_metrics[f'{name}_max'] = var_masked.max(dim=['nCells', 'nVertLevels'], skipna=True).transpose('Time',
                                                                                                               'region')
            file_metrics[f'{name}_min'] = var_masked.min(dim=['nCells', 'nVertLevels'], skipna=True).transpose('Time',
                                                                                                               'region')

            # 4. Calculate Volume-Weighted Spatial Mean
            vol_broadcasted = cell_volume_3d.expand_dims(region=iam_da.region)
            vol_masked = vol_broadcasted.where(iam_da).where(var_masked.notnull())

            weighted_sum = (var_masked * vol_masked).sum(dim=['nCells', 'nVertLevels'], skipna=True)
            total_volume = vol_masked.sum(dim=['nCells', 'nVertLevels'], skipna=True)

            file_metrics[f'{name}_mean'] = (weighted_sum / total_volume).transpose('Time', 'region')

        # Convert the dictionary of metrics into a single dataset and append to accumulator
        monthly_ds = xr.Dataset(file_metrics)
        monthly_tseries_list.append(monthly_ds)

print("Finished processing individual files. Concatenating time series...")

# Step 3: Concatenate all the processed monthly intervals along the 'Time' dimension
out_ds = xr.concat(monthly_tseries_list, dim='Time')

# Ensure correct coordinate order
out_ds = out_ds.transpose('Time', 'region')

# Save cleaner dataset to NetCDF
output_filename = f'bulk_tseries_1by1_{dx}_{sec}{subsec}.nc'
out_ds.to_netcdf(output_filename)
print(f"Successfully saved consolidated time series to: {output_filename}")