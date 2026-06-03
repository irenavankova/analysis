#!/usr/bin/env python3
import xarray as xr

# Load the dataset
ds = xr.open_dataset('/Users/ivankova/Desktop/Fris_hr/Fris_derived/nc_files/obs_tseries/obs_tseries_F1_Spin1p1.nc')

# Select site index 5 for a 2D variable
var2print = 'timeMonthly_avg_activeTracerHorMixTendency_temperatureHorMixTendency'
stprint = ds[var2print].sel(site="R02")
print(var2print)
print(stprint.values)
