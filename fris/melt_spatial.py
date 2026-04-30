#!/usr/bin/env python3

import xarray as xr
import numpy as np
import gmask_reg
import os

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#import cartopy.feature as cfeature

import cmocean
import mosaic
import matplotlib.path as mpath
from matplotlib import colors


fnum = 2
opt_plot = 3
opt_save = 1
fris_loc = '/Users/irenavankova/Desktop/Fris_hr'
#mesh_file = "/Users/irenavankova/Desktop/Fris_ncfiles/20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC08to60E3r1.spinY6_scr5.chicoma-cpu.mpaso.rst.0004-01-01_00000.nc"
#mesh_file = f'/Users/irenavankova/Desktop/Fris_ncfiles/20240227.GMPAS-JRA1p5-DIB-PISMF.TL319_FRISwISC08to60E3r1.spinY6_scr5.chicoma-cpu.mpaso.rst.0004-01-01_00000.nc'

mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}mesh.nc'
out_file = f'{fris_loc}/Fris_ncfiles/F{fnum}melt_annual_2D_Y4.nc'

iceshelves = ["Fris"]
iam = gmask_reg.get_mask(iceshelves, mesh_file)
iis = iam[0,:]

dsMesh = xr.open_dataset(mesh_file)
# Subset the dataset to ONLY include the indices in 'iam'
#dsMesh = dsMesh.isel(nCells=iis)

dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell','maxLevelCell','nEdges','nVertices','lonEdge','latEdge','lonVertex','latVertex','cellsOnEdge','cellsOnVertex','verticesOnEdge','verticesOnCell','edgesOnVertex']]
dsMesh.load()

FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
areaCell = dsMesh['areaCell'].values
lat = dsMesh['latCell'].values
lon = dsMesh['lonCell'].values
lat = lat*180/np.pi
lon = lon*180/np.pi

ds = xr.open_dataset(out_file)
#ds = ds.isel(nCells=iis)
ds = ds[['timeMonthly_avg_landIceFreshwaterFluxTotal']]

# Melt rate
rho_fw = 1000.
# seconds per day
sec_per_day = 86400.
# seconds in a year
sec_per_year = sec_per_day * 365.
lifw = ds['timeMonthly_avg_landIceFreshwaterFluxTotal']
lifw = lifw * FloatingMask * sec_per_year / rho_fw

#PLOT
if opt_plot == 1:
    plot_data = lifw.isel(Time=0).values

    plt.figure(figsize=(10, 8))

    # 2. Create the scatter plot
    # s=1 sets a small point size; cmap='RdBu_r' is good for melt/freeze
    v_min = -2
    v_max = 2

    sc = plt.scatter(lon, lat, c=plot_data, s=3, vmin=v_min, vmax=v_max, cmap='RdBu_r')

    # 3. Add a colorbar and labels
    plt.colorbar(sc, label='Freshwater Flux (m/yr)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Freshwater Flux: {iceshelves[0]}')

    # Optional: Adjust the axes to focus on the ice shelf
    plt.xlim(lon.min(), lon.max())
    plt.ylim(lat.min(), lat.max())

    plt.tight_layout()
    plt.show()

elif opt_plot == 2:
    # 1. Flatten the data and remove NaNs (important for histograms)
    # We use .values to get the numpy array from the xarray object
    data_to_plot = lifw.values.flatten()
    data_to_plot = data_to_plot[~np.isnan(data_to_plot)]

    plt.figure(figsize=(8, 6))

    # 2. Create the normalized histogram
    # density=True makes it a probability density (normalized)
    # bins=50 provides a good balance of detail
    counts, bins, patches = plt.hist(data_to_plot, bins=50, range=(-5, 5), density=True,
             color='skyblue', edgecolor='black', alpha=0.7)

    # 3. Formatting
    plt.title(f'Normalized Distribution of Freshwater Flux: {iceshelves[0]}', fontsize=12)
    plt.xlabel('Freshwater Flux (m/yr)')
    plt.ylabel('Probability Density')

    # Optional: Add a vertical line at 0 to separate melting from freezing
    plt.axvline(0, color='red', linestyle='--', linewidth=1, label='Zero Melt')
    plt.legend()

    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()


elif opt_plot == 3:
    plot_data = lifw.isel(Time=0).values

    # download and read the mesh from lcrc
    #if not os.path.exists(mesh_file):
    #    print(f"❌ Error: File not found at {mesh_file}")
    #else:
    #    dm = mosaic.datasets.open_dataset(mesh_file)
#
    #dm = dm.isel(nCells=iis)

    # define a map projection for our figure
    #projection = ccrs.SouthPolarStereo()
    projection = ccrs.PlateCarree()

    # define the transform that describes our dataset
    transform = ccrs.SouthPolarStereo()

    # create the figure and a GeoAxis
    fig, ax = plt.subplots(
        figsize=(10, 7),
        constrained_layout=True,
        subplot_kw={"projection": projection},
    )

    # create a `Descriptor` object which takes the mesh information and creates
    # the polygon coordinate arrays needed for `matplotlib.collections.PolyCollection`.
    #descriptor = mosaic.Descriptor(dsMesh, projection, transform, use_latlon=False)
    descriptor = mosaic.Descriptor(dsMesh, projection=projection)

    # using the `Descriptor` object we just created, make a pseudocolor plot of
    # the surface speed, which is defined at cell centers.
    collection = mosaic.polypcolor(
        ax,
        descriptor,
        plot_data,
        #antialiaseds=False,
        norm=colors.Normalize(vmin=-3, vmax=3),
        cmap="RdBu_r"
    )

    # Because this is not a global mesh, it's necessary to explicitly set it's extent.
    ax.set_extent([-85, -25, -84, -74], ccrs.PlateCarree())
    #ax.set_extent([-180, 180, -85, -65], ccrs.PlateCarree())

    # Below is needed for a circular boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    #ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')
    #ax.coastlines()
    fig.colorbar(collection, fraction=0.1, shrink=0.5, label="Melt rate [m/yr]")
    if opt_save == 1:
        plt.savefig(f'{fris_loc}/Fris_plots/melt/Melt_map_F{fnum}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()

