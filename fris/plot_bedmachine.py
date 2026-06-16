#!/usr/bin/env python3

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pyproj import Transformer

# -------------------------------------------------------------------------
# 1. Configuration & Setup
# -------------------------------------------------------------------------
transect_name = 'shelfbreak'  # Change this for each new track

fnameB = '/Users/ivankova/Desktop/Fris_hr/Fris_ncfiles/BedMachineAntarctica-v3.nc'
save_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/pts_4_analysis/'
n = 5  # Downsampling factor
mkm = 1000.0  # Meters to kilometers conversion
vmin = -1500

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Define the expected output file path
out_path = os.path.join(save_dir, f"{transect_name}.nc")
file_exists = os.path.exists(out_path)

# Coordinate transformer (only needed if we are generating new points)
if not file_exists:
    crs_polar = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    crs_wgs84 = "EPSG:4326"
    transformer = Transformer.from_crs(crs_polar, crs_wgs84, always_xy=True)

# -------------------------------------------------------------------------
# 2. Load BedMachine & Subset to FRIS
# -------------------------------------------------------------------------
print("Loading BedMachine...")
ds = xr.open_dataset(fnameB)

ds_fris = ds.sel(
    x=slice(-2.3e6, -0),
    y=slice(20e5, 1e5)
)

# Apply the downsampling factor to the subsetted region
ds_sub = ds_fris.isel(x=slice(0, None, n), y=slice(0, None, n))

# Convert subset axes coordinates to kilometers for plotting
x_km = ds_sub['x'] / mkm
y_km = ds_sub['y'] / mkm

print("Loading values...")
var2plot = ds_sub['bed'].values

# -------------------------------------------------------------------------
# 3. Setup Plotting Environment
# -------------------------------------------------------------------------
if not file_exists:
    plt.ion()  # Interactive mode ON if we need to click points
else:
    plt.ioff()  # Interactive mode OFF if we are just loading

fig, ax = plt.subplots(figsize=(14, 8))
X_grid, Y_grid = np.meshgrid(x_km, y_km)

# Plot background bed data
pcm = ax.pcolormesh(X_grid, Y_grid, var2plot, cmap='pink', vmin=vmin, vmax=0, shading='nearest')
fig.colorbar(pcm, ax=ax, label='Water Column Thickness [m]')

ax.set_xlabel("X coordinate [km]")
ax.set_ylabel("Y coordinate [km]")

# -------------------------------------------------------------------------
# 4. Conditional Logic: Load File OR Interactive Selection
# -------------------------------------------------------------------------
if file_exists:
    # --- ROUTINE A: LOAD EXISTING FILE ---
    print(f"\n[FOUND] Existing file detected: {out_path}")
    print("Loading coordinates from file...")

    ds_points = xr.open_dataset(out_path)
    click_x_km = ds_points['x'].values / mkm
    click_y_km = ds_points['y'].values / mkm

    ax.set_title(f"FRIS Region - {transect_name}\nLoaded from {transect_name}.nc")
    ax.plot(click_x_km, click_y_km, 'g-o', linewidth=2, label=f"{transect_name} (Loaded)")

else:
    # --- ROUTINE B: INTERACTIVE SELECTION (ginput) ---
    print(f"\n[NOT FOUND] No file at: {out_path}")
    print("--> CLICK ON THE MAP NOW. Press ENTER when done.")

    ax.set_title(f"FRIS Region - {transect_name}\nLEFT CLICK to add points. Press ENTER/RETURN when finished.")
    plt.draw()

    clicked_points = plt.ginput(n=-1, timeout=0)

    if not clicked_points:
        print("No points selected. Exiting script.")
        plt.close(fig)
        exit()

    clicked_points = np.array(clicked_points)
    click_x_km = clicked_points[:, 0]
    click_y_km = clicked_points[:, 1]

    ax.plot(click_x_km, click_y_km, 'r-o', linewidth=2, label=f"{transect_name} (New)")

    # Process coordinates and save to NetCDF (Only if newly created)
    x_meters = click_x_km * mkm
    y_meters = click_y_km * mkm
    longitude, latitude = transformer.transform(x_meters, y_meters)
    nPoints = len(latitude)

    dataset_out = xr.Dataset(
        data_vars=dict(
            x=(["nPoints"], x_meters),
            y=(["nPoints"], y_meters),
            longitude=(["nPoints"], longitude),
            latitude=(["nPoints"], latitude)
        ),
        coords=dict(nPoints=np.arange(nPoints)),
        attrs=dict(
            description='Manually drawn FRIS transect via Python ginput',
            transect_name=transect_name
        ),
    )
    dataset_out.to_netcdf(out_path)
    print(f"Successfully saved coordinates to: {out_path}")

# -------------------------------------------------------------------------
# NEW: Add Site Labels to the Plot
# -------------------------------------------------------------------------
# Adjust the offset values to shift the text safely away from the dot
xy_offset = 15  # km distance

for i, (x, y) in enumerate(zip(click_x_km, click_y_km)):
    ax.text(
        x + xy_offset,
        y + xy_offset,
        f"{i+1}",
        fontsize=8,
        fontweight='bold',
        color='white',
        bbox=dict(facecolor='none', alpha=0.6, edgecolor='none', pad=1) # Makes labels readable over dark/pink spots
    )

# -------------------------------------------------------------------------
# 5. Finalize and Save Plot Image
# -------------------------------------------------------------------------
ax.legend()
plt.draw()

plot_save_path = os.path.join(save_dir, f"{transect_name}_plot.png")
fig.savefig(plot_save_path, bbox_inches='tight', dpi=300)
print(f"Successfully saved plot to: {plot_save_path}")

# Keep plot window open for review
if not file_exists:
    plt.ioff()
plt.show()