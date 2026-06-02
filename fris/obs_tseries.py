#!/usr/bin/env python3

import glob
import xarray as xr
import numpy as np

# Import the coordinates directly from your existing script
from fris_coordinates import sites_config, find_nearest_mpas_cells

# --- Run Configuration ---
Fnum = '4'
dx = f'F{Fnum}'
sec = 'Spin6'
subsec = 'p1'

run_name = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"
fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'
fpath_mask = fpath
run_name_mask = run_name
mask_file = f'{fpath_mask}/{run_name_mask}.mpaso.rst.0002-01-01_00000.nc'

tot_avail = 0
if sec == 'Spin6':
    subsec = ''
    tot_avail = 1

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

with xr.open_dataset(mask_file) as dsMesh:
    site_da = find_nearest_mpas_cells(sites_config, dsMesh)

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
    'timeMonthly_avg_potentialDensity'
]

vars_2d = [
    'timeMonthly_avg_landIceFreshwaterFlux',
    'timeMonthly_avg_landIceFreshwaterFluxTotal',
    'timeMonthly_avg_ssh'
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

        # Pull 3D layered fields safely (Vertical dimension nVertLevels is fully untouched)
        for varname in vars_3d:
            if varname in ds:
                file_metrics[varname] = ds[varname].isel(nCells=site_da)

        # Bundle file metrics together, preserving native time frames
        monthly_ds = xr.Dataset(file_metrics)
        monthly_tseries_list.append(monthly_ds)

print("Concatenating monthly site profiles into master dataset...")
out_ds = xr.concat(monthly_tseries_list, dim='Time')

# Ensure tidy dimension tracking order (Time, site, nVertLevels)
out_ds = out_ds.transpose('Time', 'site', ...)

# Save cleaner dataset to NetCDF
output_filename = f'site_tseries_{dx}_{sec}{subsec}.nc'
out_ds.to_netcdf(output_filename)
print(f"Successfully consolidated point profiles to: {output_filename}")