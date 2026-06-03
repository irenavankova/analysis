#!/usr/bin/env python3

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.signal import medfilt2d  # Added for medfilt2 equivalent

# -------------------------------------------------------------------------
# 1. Configuration & Setup
# -------------------------------------------------------------------------
transect_name = 'ronnecenter'  # Change this for each new track

fnameB = '/Users/ivankova/Desktop/Fris_hr/Fris_ncfiles/BedMachineAntarctica-v3.nc'
save_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/pts_4_analysis/'
n = 5                                # Downsampling factor (reduced since region is smaller)
mkm = 1000.0                         # Meters to kilometers conversion
vmin = -1500

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Coordinate transformer for Polar Stereo (m) -> Lon/Lat (degrees)
crs_polar = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_wgs84 = "EPSG:4326"
transformer = Transformer.from_crs(crs_polar, crs_wgs84, always_xy=True)

# -------------------------------------------------------------------------
# 2. Load BedMachine & Subset to FRIS
# -------------------------------------------------------------------------
print("Loading BedMachine...")
ds = xr.open_dataset(fnameB)

# BedMachine Antarctica Y coordinates go from positive to negative (descending).
# Therefore, slice must go from LARGER to SMALLER number.
ds_fris = ds.sel(
    #x=slice(-1.6e6, -0.5e6),
    #y=slice(11e5, 1e5)        # Swapped: 11e5 first, then 1e5
    x=slice(-2.3e6, -0),
    y=slice(20e5, 1e5)        # Swapped: 11e5 first, then 1e5
)

# Apply the downsampling factor to the subsetted region
ds_sub = ds_fris.isel(x=slice(0, None, n), y=slice(0, None, n))

# Convert subset axes coordinates to kilometers for plotting
x_km = ds_sub['x'] / mkm
y_km = ds_sub['y'] / mkm

print("Loading values...")
# Extract required fields as numpy arrays
#surface = ds_sub['surface'].values
#thickness = ds_sub['thickness'].values
bed = ds_sub['bed'].values
#mask_sub = ds_sub['mask'].values

# Calculate raw water column thickness
#wct_raw = (surface - thickness) - bed

print("Loading values...")
var2plot = bed

# -------------------------------------------------------------------------
# 3. Plot & Interactive Selection (ginput)
# -------------------------------------------------------------------------
plt.ion()
fig, ax = plt.subplots(figsize=(11, 10))

print("Plotting FRIS regional map...")

# Create 2D coordinate grids matching the shape of wct_sub
X_grid, Y_grid = np.meshgrid(x_km, y_km)

# Plot using the 2D grids (use shading='nearest' to avoid boundary size issues)
pcm = ax.pcolormesh(X_grid, Y_grid, var2plot, cmap='pink', vmin=vmin, vmax=0, shading='nearest')
fig.colorbar(pcm, ax=ax, label='Water Column Thickness [m]')

# Update the contour line to use the 2D grid as well
#ax.contour(X_grid, Y_grid, mask_sub, levels=[1, 2], colors='w', linewidths=0.75)

ax.set_title(f"FRIS Region - {transect_name}\nLEFT CLICK to add points. Press ENTER/RETURN when finished.")
ax.set_xlabel("X coordinate [km]")
ax.set_ylabel("Y coordinate [km]")
plt.draw()

# Call ginput
print("\n--> CLICK ON THE MAP NOW. Press ENTER when done.")
clicked_points = plt.ginput(n=-1, timeout=0)

if not clicked_points:
    print("No points selected. Exiting script.")
    plt.close(fig)
    exit()

# Extract clicked x and y (currently in km)
clicked_points = np.array(clicked_points)
click_x_km = clicked_points[:, 0]
click_y_km = clicked_points[:, 1]

# Show your track on the map as a confirmation visual
ax.plot(click_x_km, click_y_km, 'r-o', linewidth=2, label=transect_name)
ax.legend()
plt.draw()
print(f"Captured {len(clicked_points)} points.")

# --- ADDED: Save the plot to the same directory as the input .nc file ---
plot_save_path = os.path.join(save_dir, f"{transect_name}_plot.png")
fig.savefig(plot_save_path, bbox_inches='tight', dpi=300)
print(f"Successfully saved plot to: {plot_save_path}")
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# 4. Process Coordinates & Save Directly to NetCDF
# -------------------------------------------------------------------------
# Convert the clicked km coordinates back to meters for pyproj
x_meters = click_x_km * mkm
y_meters = click_y_km * mkm

# Transform polar stereographic meters back to standard Lat/Lon
longitude, latitude = transformer.transform(x_meters, y_meters)

nPoints = len(latitude)

# Bundle everything cleanly inside an xarray Dataset
dataset_out = xr.Dataset(
    data_vars=dict(
        x=(["nPoints"], x_meters),
        y=(["nPoints"], y_meters),
        longitude=(["nPoints"], longitude),
        latitude=(["nPoints"], latitude)
    ),
    coords=dict(
        nPoints=np.arange(nPoints)
    ),
    attrs=dict(
        description='Manually drawn FRIS transect via Python ginput',
        transect_name=transect_name
    ),
)

# Save the final NetCDF
out_path = os.path.join(save_dir, f"{transect_name}.nc")
dataset_out.to_netcdf(out_path)
print(f"\nSuccessfully saved coordinates to: {out_path}")

# Keep plot open for a moment to review before script terminates
plt.ioff()
plt.show()