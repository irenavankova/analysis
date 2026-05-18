#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

Fnum = '8'
dx = f'F{Fnum}'
sec = 'Spin6'

run_name = f"20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC0{Fnum}to60E3r1.spinY6_scr5.chicoma-cpu"

fpath = f'/pscratch/sd/v/vankova/lanl/FRIS_Irena/FRIS_spinY6/{run_name}/run'

# Step 1: Load the mask and grid metrics from the restart NetCDF file
mask_file = f'{fpath}/{run_name}.mpaso.rst.0002-01-01_00000.nc'

mask_ds = xr.open_dataset(mask_file)
areaCell = mask_ds['areaCell']  # Horizontal grid cell area (nCells,)

regs = ["FRIS", "FRISshelf", "RonneDshelf", "FilchnerDshelf", "BerknerBank", "RonneDcavity", "FilchnerDcavity", "BerknerSouth"]
iam = gmask_reg.get_mask(regs, mask_file)


# Step 2: Update preprocess to load only the necessary 3D tracer variables
def preprocess(dss):
    return dss[
        [
            "timeMonthly_avg_activeTracers_temperature",
            "timeMonthly_avg_activeTracers_salinity",
            "timeMonthly_avg_potentialDensity",
            "timeMonthly_avg_layerThickness"  # Essential for 3D cell volume weighting
        ]
    ]


ds = xr.open_mfdataset(f"{fpath}/{run_name}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc",
                           combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True,
                           preprocess=preprocess)

# --- Prepare 3D Variables and Volumetric Weights ---
T_3d = ds['timeMonthly_avg_activeTracers_temperature']
S_3d = ds['timeMonthly_avg_activeTracers_salinity']
Rho_3d = ds['timeMonthly_avg_potentialDensity']

# Calculate 3D Cell Volume: areaCell (nCells) * layerThickness (Time, nCells, nVertLevels)
h_3d = ds['timeMonthly_avg_layerThickness']
cell_volume_3d = areaCell * h_3d

# --- Vectorized Masking Setup ---
# Reshape the region mask into an xarray DataArray for seamless broadcasting: (region, nCells)
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': regs})

# Dictionary to collect all processed time series arrays
tseries_dict = {}

# --- 3D Quantities Masking and Reductions ---
vars_3d = {
    'temp': T_3d,
    'salt': S_3d,
    'rho': Rho_3d
}

for name, var in vars_3d.items():
    # 1. Broadcast the 3D variable across the horizontal region mask dimensions
    # New shape: (region, Time, nCells, nVertLevels)
    var_broadcasted = var.expand_dims(region=iam_da.region)

    # 2. Mask out data points falling outside the horizontal boundaries of each region
    # Additionally mask out exact 0.0 or NaN fills common in deeper uninitialized layers
    var_masked = var_broadcasted.where(iam_da).where(var_broadcasted.notnull())

    # 3. Calculate spatial Minimum and Maximum over both spatial dimensions (nCells, nVertLevels)
    # skipna=True guarantees that NaNs do not propagate and cause the output to be NaN
    tseries_dict[f'{name}_max'] = var_masked.max(dim=['nCells', 'nVertLevels'], skipna=True).transpose('Time', 'region')
    tseries_dict[f'{name}_min'] = var_masked.min(dim=['nCells', 'nVertLevels'], skipna=True).transpose('Time', 'region')

    # 4. Calculate Volume-Weighted Spatial Mean
    # Ensure the cell volume weights mirror the region masking and handle missing values identically
    vol_broadcasted = cell_volume_3d.expand_dims(region=iam_da.region)
    vol_masked = vol_broadcasted.where(iam_da).where(var_masked.notnull())

    # skipna=True handles unexpected internal NaNs elegantly during summation steps
    weighted_sum = (var_masked * vol_masked).sum(dim=['nCells', 'nVertLevels'], skipna=True)
    total_volume = vol_masked.sum(dim=['nCells', 'nVertLevels'], skipna=True)

    # Final Division (invalid cells without volume evaluate to NaN naturally here)
    tseries_dict[f'{name}_mean'] = (weighted_sum / total_volume).transpose('Time', 'region')

# --- Save cleaner dataset to NetCDF ---
out_ds = xr.Dataset(tseries_dict)
out_ds.to_netcdf(f'ocean_reg_tseries_{dx}_{sec}.nc')