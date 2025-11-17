#!/usr/bin/env python3

import xarray as xr
import numpy as np

import gmask_reg

ymem = '701'

fpath_pismf = '/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_ensmean/run/'

if ymem == 'ave':
    fpath = f'{fpath_pismf}'
else:
    fpath = f'/lcrc/group/acme/ac.dcomeau/scratch/chrys/E3SMv2_1/v2_1.SORRM.ssp370_0{ymem}/'
    fpath = f'{fpath}archive/ocn/hist/'

secPerYear = 365 * 24 * 60 * 60
kgingt = 1e12

# Step 1: Load the mask from a different NetCDF file
mask_file = f'{fpath_pismf}v2_1.SORRM.ssp370_ensmean.mpaso.rst.2015-01-01_00000.nc'  # Specify the mask file
mask_ds = xr.open_dataset(mask_file)
landIceFloatingMask = mask_ds['landIceFloatingMask']  # Assuming the variable name is 'mask' and it's of shape (ncells,)
areaCell = mask_ds['areaCell']  # Assuming the variable name is 'mask' and it's of shape (ncells,)

iceshelves = ["Antarctica", "Belli", "Amundsen", "Ross", "Eastant", "Amery", "Dml", "Fris", "Larsens"]
iam = gmask_reg.get_mask(iceshelves, mask_file)

# Step 2: Open the 100 NetCDF files and concatenate them along the 'time' dimension
def preprocess(dss):
    return dss[
        [
            "timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature",
            "timeMonthly_avg_velocityZonalSquared",
            "timeMonthly_avg_velocityMeridionalSquared",
            "timeMonthly_avg_activeTracers_temperature"
        ]
    ]

if ymem == 'ave':
    ds = xr.open_mfdataset(f"{fpath}v2_1.SORRM.ssp370_ensmean.mpaso.hist.am.timeSeriesStatsMonthly.*.nc",
                           combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True,
                           preprocess=preprocess)
else:
    ds = xr.open_mfdataset(f"{fpath}v2_1.SORRM.ssp370_0{ymem}.mpaso.hist.am.timeSeriesStatsMonthly.*.nc", combine='by_coords', chunks={'time': 12}, parallel=True, decode_timedelta=True, preprocess=preprocess)


landIceFloatingMask = landIceFloatingMask.squeeze('Time')

# Melt rate
#lifw = ds['timeMonthly_avg_landIceFreshwaterFlux']
#lifw = lifw * landIceFloatingMask * areaCell * secPerYear / kgingt

Ti = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']
Vsq = ds['timeMonthly_avg_velocityMeridionalSquared'].isel(nVertLevels=0)
Usq = ds['timeMonthly_avg_velocityZonalSquared'].isel(nVertLevels=0)
T = ds['timeMonthly_avg_activeTracers_temperature'].isel(nVertLevels=0)

# Access the global attribute
Cd = ds.attrs['config_land_ice_flux_topDragCoeff']
utidal = ds.attrs['config_land_ice_flux_rms_tidal_velocity']

utsq = np.sqrt(Cd * (Usq + Vsq + utidal**2)) * (T - Ti)

utsq = utsq * landIceFloatingMask * areaCell * secPerYear / kgingt

# --- Vectorized masking over all regions ---
# Stack region mask into xarray DataArray for better alignment
iam_da = xr.DataArray(iam, dims=['region', 'nCells'], coords={'region': iceshelves})

# Expand lifw to shape (region, Time, nCells) by broadcasting
utsq_broadcasted = utsq.expand_dims(region=iam_da.region)

# Use masking across all regions simultaneously
utsq_masked = utsq_broadcasted.where(iam_da)

# Sum over nCells → result shape: (region, Time)
utsq_tseries = utsq_masked.sum(dim='nCells')

# Transpose to have (Time, region) if desired
utsq_tseries = utsq_tseries.transpose('Time', 'region')
utsq_tseries.name = 'utsq'

# Save to NetCDF
utsq_tseries.to_dataset().to_netcdf(f'utsq_reg_pismf_{ymem}_tseries.nc')

