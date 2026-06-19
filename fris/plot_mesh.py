#!/usr/bin/env python3

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mosaic
from matplotlib import colors

# --- Configurations ---
opt_save = 0
fris_loc = '/Users/ivankova/Desktop/Fris_hr'
fnum = 8
mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'

# --- Load Mesh Dataset ---
dsMesh = xr.open_dataset(mesh_file)
dsMesh = dsMesh[
    ['latCell', 'lonCell', 'areaCell', 'nEdges', 'nVertices', 'lonEdge',
     'latEdge', 'lonVertex', 'latVertex', 'cellsOnEdge', 'cellsOnVertex',
     'verticesOnEdge', 'verticesOnCell', 'edgesOnVertex']
]
dsMesh.load()

# --- Extract & Process Target Variable ---
# Convert area from square meters to square kilometers, then take the square root
plot_data = np.sqrt(dsMesh['areaCell'].values / 1e6)
plotlabel = "Cell Area [km$^2$]"

# Normalize dynamically based on the grid resolution range
v_min = 0.5
v_max = 13.0
norm = colors.Normalize(vmin=v_min, vmax=v_max)
varcmap = "gist_earth"  # Good sequential colormap for grid cell sizes
# Replace the old norm with a BoundaryNorm for exact interval matching
specified_ticks = [1, 2, 4, 8, 12]
#norm = colors.BoundaryNorm(boundaries=specified_ticks, ncolors=plt.get_cmap(varcmap).N)

# --- Setup Plot ---
projection = ccrs.SouthPolarStereo(central_longitude=-75)
fig, ax = plt.subplots(
    figsize=(5, 4.5),
    constrained_layout=True,
    subplot_kw={"projection": projection},
)

# Render MPAS unstructured mesh polygons
descriptor = mosaic.Descriptor(dsMesh, projection=projection)
collection = mosaic.polypcolor(
    ax,
    descriptor,
    plot_data,
    norm=norm,
    cmap=varcmap,
    linewidth = 0.0,
    edgecolors='face'
    #antialiased=False
)
#mosaic.coastlines(ax, descriptor)
# --- Map Formatting ---
ax.set_extent([-80, 0, -84, -64], ccrs.PlateCarree()) #GOOD
#ax.set_extent([-82, -25, -81, -72], ccrs.PlateCarree())
#ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
ax.set_aspect('auto')
ax.gridlines(draw_labels=True)
ax.set_facecolor('lightgray')

# --- Add Colorbar ---
# Create the colorbar and enforce explicit limits on its axes
cb = fig.colorbar(collection, fraction=0.1, shrink=1.0, label=plotlabel, extend='both',ticks=specified_ticks)
cb.ax.set_ylim(v_min, v_max)

# --- Save or Show ---
if opt_save == 1:
    plt.savefig(f'{fris_loc}/Fris_plots/mesh/sqrtAreaCell_F{fnum}.png', bbox_inches='tight', dpi=300)
else:
    plt.show()