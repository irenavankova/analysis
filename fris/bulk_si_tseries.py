#!/usr/bin/env python3

import glob
import xarray as xr
import numpy as np
import gmask_reg

Fnum = '8'
dx = f'F{Fnum}'
sec = 'Spin6'
subsec = 'p1'
opt_noGL = 1

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

# Scale factor to convert m/s to m/yr
units_scale_factor = 60 * 60 * 24 * 365

# Step 1: Load the mask and grid metrics from the ocean restart file
mask_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

mask_ds = xr.open_dataset(mask_file)
areaCell = mask_ds['areaCell']  # Horizontal grid cell area (nCells,)

regs = ["FRISshelf", "RonneDshelf", "FilchnerDshelf", "BerknerBank"]

print(f"Getting masks...")
iam = gmask_reg.get_mask(regs, mask_file, opt_noGL=opt_noGL)

print(f"Reshaping masks...")
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': regs})

print(f"Getting a file list...")
file_pattern = f"{fpath}/{run_name}.mpassi.hist.am.timeSeriesStatsMonthly.*.nc"
file_list = sorted(glob.glob(file_pattern))

if not file_list:
    raise FileNotFoundError(f"No history files found matching pattern: {file_pattern}")

# --- NEW STEP: Pre-screening files to find and exclude cells where sea-ice volume > 100 ---
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
# ------------------------------------------------------------------------------------------

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

        # 2. Extract and Calculate Total Sea Ice Production
        congelation = ds['timeMonthly_avg_congelation']
        frazil = ds['timeMonthly_avg_frazilFormation']
        snowice = ds['timeMonthly_avg_snowiceFormation']
        ice_production = (congelation + frazil + snowice) * units_scale_factor

        # 3. Extract and Calculate Total Sea Ice Melting
        basal = ds['timeMonthly_avg_basalIceMelt']
        surface = ds['timeMonthly_avg_surfaceIceMelt']
        lateral = ds['timeMonthly_avg_lateralIceMelt']
        ice_melting = (basal + surface + lateral) * units_scale_factor

        file_metrics = {}

        # Dictionary of fields where we want area-weighted means
        vars_for_mean = {
            'ice_concentration': ice_concentration,
            'ice_thickness': ice_thickness,
            'ice_production': ice_production,
            'ice_melting': ice_melting
        }

        # Broadcast grid area safely across regions and apply masks once per file loop
        area_broadcasted = areaCell.expand_dims(region=iam_da.region)
        area_masked = area_broadcasted.where(iam_da).where(valid_cells_mask)

        # ---- Calculate Area-Weighted Means ----
        for name, var in vars_for_mean.items():
            var_broadcasted = var.expand_dims(region=iam_da.region)
            var_masked = var_broadcasted.where(iam_da).where(valid_cells_mask)

            # Ensure area mask only counts cells where the data variable itself isn't NaN
            active_area = area_masked.where(var_masked.notnull())

            # Area-weighted mean: sum(var * area) / sum(area)
            weighted_sum = (var_masked * active_area).sum(dim='nCells', skipna=True)
            total_area = active_area.sum(dim='nCells', skipna=True)

            file_metrics[f'mean_{name}'] = (weighted_sum / total_area).transpose('Time', 'region')

        # ---- Calculate Integrated Ice Volume ----
        # Volume = thickness * cell_area.
        thick_broadcasted = ice_thickness.expand_dims(region=iam_da.region)
        thick_masked = thick_broadcasted.where(iam_da).where(valid_cells_mask)

        # Integrated volume over the valid grid region
        file_metrics['integrated_ice_volume'] = (thick_masked * area_masked).sum(dim='nCells', skipna=True).transpose(
            'Time', 'region')

        # ---- Calculate Maximum Ice Thickness ----
        thick_masked_allValid = thick_broadcasted.where(iam_da)
        file_metrics['max_ice_thickness'] = thick_masked_allValid.max(dim='nCells', skipna=True).transpose('Time', 'region')

        # Convert the dictionary of metrics into a single dataset and append to accumulator
        monthly_ds = xr.Dataset(file_metrics)
        monthly_tseries_list.append(monthly_ds)

print("Finished processing individual files. Concatenating time series...")

# Step 3: Concatenate all the processed monthly intervals along the 'Time' dimension
out_ds = xr.concat(monthly_tseries_list, dim='Time')

# Ensure correct coordinate order
out_ds = out_ds.transpose('Time', 'region')

# Save cleaner dataset to NetCDF
output_filename = f'bulk_seaice_tseries_{dx}_{sec}{subsec}.nc'
out_ds.to_netcdf(output_filename)
print(f"Successfully saved consolidated time series to: {output_filename}")