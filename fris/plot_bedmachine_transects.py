#!/usr/bin/env python3

import os
import json
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pyproj import Transformer
from matplotlib.lines import Line2D

# -------------------------------------------------------------------------
# 1. Configuration & Setup
# -------------------------------------------------------------------------
fnameB = '/Users/ivankova/Desktop/Fris_hr/Fris_ncfiles/BedMachineAntarctica-v3.nc'
base_save_dir = '/Users/ivankova/Desktop/Fris_hr/Fris_derived/pts_transects/'

n = 5  # Downsampling factor
mkm = 1000.0  # Meters to kilometers conversion
vmin = 0
vmax = 1500

crs_polar = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_wgs84 = "EPSG:4326"
transformer = Transformer.from_crs(crs_polar, crs_wgs84, always_xy=True)

# -------------------------------------------------------------------------
# 2. Load BedMachine & Subset to FRIS (Optimized)
# -------------------------------------------------------------------------
print("Loading BedMachine...")
ds = xr.open_dataset(fnameB)

# Subset to the FRIS region and immediately downsample
ds_sub = ds.sel(
    x=slice(-2.3e6, -0),
    y=slice(20e5, 1e5)
).isel(
    x=slice(0, None, n),
    y=slice(0, None, n)
)

print("Loading localized values efficiently...")
# Calculate water thickness lazily inside xarray before converting to numpy
water_column_thickness = ds_sub['surface'] - ds_sub['thickness'] - ds_sub['bed']
var2plot = water_column_thickness.values

# Convert subset axes coordinates to kilometers for plotting
x_km = ds_sub['x'].values / mkm
y_km = ds_sub['y'].values / mkm

# -------------------------------------------------------------------------
# 3. Setup Plotting Environment & Interactive Selection
# -------------------------------------------------------------------------
plt.ion()  # Interactive mode ON
fig, ax = plt.subplots(figsize=(14, 8))
X_grid, Y_grid = np.meshgrid(x_km, y_km)

# Plot background bed data
pcm = ax.pcolormesh(X_grid, Y_grid, var2plot, cmap='pink', vmin=vmin, vmax=vmax, shading='nearest')
fig.colorbar(pcm, ax=ax, label='Water Column Thickness [m]')

ax.set_xlabel("X coordinate [km]")
ax.set_ylabel("Y coordinate [km]")

print("\n--> CLICK AN EVEN NUMBER OF POINTS ON THE MAP NOW. Press ENTER when finished.")
ax.set_title(
    "FRIS Region - Multi-Transect Input\nLEFT CLICK an EVEN number of points (Pairs). Press ENTER/RETURN when finished.")
plt.draw()

clicked_points = plt.ginput(n=-1, timeout=0)

# CRITICAL CHECK: Enforce an even number of points greater than zero
num_points = len(clicked_points)
if num_points == 0 or num_points % 2 != 0:
    plt.close(fig)
    raise ValueError(
        f"Error: You must select an EVEN number of points. You provided {num_points}. No files were saved.")

clicked_points = np.array(clicked_points)
total_transects = num_points // 2
print(f"Detected {num_points} points. Processing {total_transects} transects...")

# -------------------------------------------------------------------------
# 4. Loop Through Pairs to Process & Save Individual Transects
# -------------------------------------------------------------------------
xy_offset = 15  # km label text offset

for t_idx in range(total_transects):
    # Calculate zero-padded transect ID (F01, F02, ...)
    transect_num = t_idx + 1
    transect_name = f"F{transect_num:02d}"

    # Establish folder path for this specific transect
    transect_save_dir = os.path.join(base_save_dir, transect_name)
    if not os.path.exists(transect_save_dir):
        os.makedirs(transect_save_dir)

    out_path = os.path.join(transect_save_dir, f"{transect_name}.nc")
    geojson_path = os.path.join(transect_save_dir, "transect.geojson")

    # Extract the pair of points for this transect
    start_idx = t_idx * 2
    end_idx = start_idx + 2

    click_x_km = clicked_points[start_idx:end_idx, 0]
    click_y_km = clicked_points[start_idx:end_idx, 1]

    # --- Plotting formatting for this segment ---
    # Draw line connecting them in black
    ax.plot(click_x_km, click_y_km, color='black', linewidth=2, zorder=2, label='_nolegend_')
    # First point in red, second point in green
    ax.scatter(click_x_km[0], click_y_km[0], color='red', s=50, zorder=3, label='_nolegend_')
    ax.scatter(click_x_km[1], click_y_km[1], color='green', s=50, zorder=3, label='_nolegend_')

    # Add point numbering labels next to the markers
    for i, (x, y) in enumerate(zip(click_x_km, click_y_km)):
        ax.text(
            x + xy_offset,
            y + xy_offset,
            f"{transect_name} P_{i + 1}",
            fontsize=8,
            fontweight='bold',
            color='white',
            bbox=dict(facecolor='none', alpha=0.6, edgecolor='none', pad=1)
        )

    # --- Data Processing and Conversions ---
    x_meters = click_x_km * mkm
    y_meters = click_y_km * mkm
    longitude, latitude = transformer.transform(x_meters, y_meters)

    # 1. Save NetCDF
    dataset_out = xr.Dataset(
        data_vars=dict(
            x=(["nPoints"], x_meters),
            y=(["nPoints"], y_meters),
            longitude=(["nPoints"], longitude),
            latitude=(["nPoints"], latitude)
        ),
        coords=dict(nPoints=np.arange(2)),
        attrs=dict(
            description=f'Manually drawn FRIS transect {transect_name} via Python multi-pair script',
            transect_name=transect_name
        ),
    )
    dataset_out.to_netcdf(out_path)
    print(f"[{transect_name}] Saved NetCDF to: {out_path}")

    # 2. Save GeoJSON
    geojson_coordinates = [[lon, lat] for lon, lat in zip(longitude, latitude)]
    geojson_data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "name": transect_name,
                    "tags": "manually_drawn_sections",
                    "object": "transect",
                    "component": "ocean",
                    "author": "Python multi-pair ginput script",
                    "note": f"FRIS region generated transect segment: {transect_name}"
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": geojson_coordinates
                }
            }
        ]
    }

    with open(geojson_path, 'w') as f:
        json.dump(geojson_data, f, indent=4)
    print(f"[{transect_name}] Saved GeoJSON to: {geojson_path}")

# -------------------------------------------------------------------------
# 5. Finalize and Save Collective Plot Image to main directory
# -------------------------------------------------------------------------
legend_elements = [
    Line2D([0], [0], color='black', lw=2, label='Transect Path'),
    Line2D([0], [0], marker='o', color='none', markerfacecolor='red', markersize=8, label='Start Point (1)'),
    Line2D([0], [0], marker='o', color='none', markerfacecolor='green', markersize=8, label='End Point (2)')
]
ax.legend(handles=legend_elements)

ax.set_title(f"FRIS Region - All Generated Transects (Total: {total_transects})")
plt.draw()

# Save the full multi-track plot to the main base directory
if not os.path.exists(base_save_dir):
    os.makedirs(base_save_dir)

plot_save_path = os.path.join(base_save_dir, "all_transects_plot.png")
fig.savefig(plot_save_path, bbox_inches='tight', dpi=300)
print(f"\n[SUCCESS] Saved collective plot figure to: {plot_save_path}")

# Keep plot window open for review
plt.ioff()
plt.show()