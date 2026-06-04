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
opt_noGL = 0
fy = '2-4'
opt_save = 1
#variables_to_plot = ['lifw', 'Tstar', 'Ti', 'Tbl', 'ustar']
#variables_to_plot = ['Ti', 'Tbl']
variables_to_plot = ['Sbot']

fris_loc = '/Users/ivankova/Desktop/Fris_hr'

mesh_file = f'/Users/ivankova/Desktop/Fris_hr/E3SM_init/ocean.ECwISC30to60E2r1.230220.nc'
out_file = mesh_file

#mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'
#out_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}melt_annual_2D_Y{fy}.nc'

opt_region = False

iceshelves = ["Shelf"]
if opt_region == True:
    iam = gmask_reg.get_mask(iceshelves, mesh_file, opt_noGL=opt_noGL, opt_wct=1)
    iis = iam[0,:]


dsMesh = xr.open_dataset(mesh_file)

dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell','maxLevelCell','nEdges','nVertices','lonEdge','latEdge','lonVertex','latVertex','cellsOnEdge','cellsOnVertex','verticesOnEdge','verticesOnCell','edgesOnVertex']]
dsMesh.load()

FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
areaCell = dsMesh['areaCell'].values
lat = dsMesh['latCell'].values
lon = dsMesh['lonCell'].values
lat = lat*180/np.pi
lon = lon*180/np.pi

max_level_cell = dsMesh['maxLevelCell'].values - 1  # Convert from 1-indexed to 0-indexed

ds = xr.open_dataset(out_file)
#ds = ds[['timeMonthly_avg_landIceFreshwaterFluxTotal']]

# Melt rate
rho_fw = 1000.
sec_per_day = 86400.
sec_per_year = sec_per_day * 365.
mtocm = 100

# Start the for loop wrapper
for var2plot in variables_to_plot:
    print(f"Processing plot for: {var2plot}")

    if var2plot == 'lifw':
        lifw = ds['timeMonthly_avg_landIceFreshwaterFluxTotal']
        lifw = lifw * FloatingMask * sec_per_year / rho_fw
        plot_data = lifw.isel(Time=0).values
        v_min = -3
        v_max = 3
        plotlabel = "Melt rate [m/yr]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    if var2plot == 'freeze':
        lifw = ds['timeMonthly_avg_landIceFreshwaterFluxTotal']
        lifw = lifw * FloatingMask * sec_per_year / rho_fw
        plot_data = lifw.isel(Time=0).values
        plot_data[plot_data > 0] = np.nan
        v_min = -1000
        v_max = 0
        plotlabel = "Freezing rate [m/yr]"
        varcmap = "pink"
        norm = colors.AsinhNorm(linear_width=1.0, vmin=-200, vmax=0)
    elif var2plot == 'Tstar':
        Tcb = ds['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature']  # Use the first time step (Time = 1)
        Tci = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']  # Use the first time step (Time = 1)
        Tc = (Tcb-Tci)*FloatingMask
        plot_data = Tc.isel(Time=0).values
        v_min = -1
        v_max = 1
        plotlabel = "Thermal driving [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Ti':
        Tci = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']  # Use the first time step (Time = 1)
        plot_data = Tci.isel(Time=0).values
        v_min = -2.8
        v_max = -1.9
        plotlabel = "Temperature interface [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Tbl':
        Tcb = ds[
            'timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature']  # Use the first time step (Time = 1)
        plot_data = Tcb.isel(Time=0).values
        v_min = -2.8
        v_max = -1.9
        plotlabel = "Temperature BL [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'ustar':
        Uc = ds['timeMonthly_avg_landIceFrictionVelocity']  # Use the first time step (Time = 1)
        #Uc = Uc*mtocm*FloatingMask
        Uc = np.where(FloatingMask == 1, Uc * mtocm, np.nan)
        #plot_data = Uc.isel(Time=0).values
        plot_data = Uc
        v_min = 0.4
        v_max = 1.2
        plotlabel = "Friction velocity [cm/s]"
        varcmap = "hot_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Sbot':
        # Assuming the variable name in your nc file is something like 'timeMonthly_avg_salinity'
        # Adjust the string key to match your actual NetCDF file variable name
        if 'timeMonthly_avg_activeTracers_salinity' in ds:
            salinity_3d = ds['timeMonthly_avg_activeTracers_salinity'].isel(Time=0).values
        elif 'salinity' in ds:
            salinity_3d = ds['salinity'].isel(Time=0).values
        else:
            raise KeyError(
                "Neither 'timeMonthly_avg_activeTracers_salinity' nor 'activeTracers_salinity' was found in the dataset.")

        # salinity_3d has shape (nCells, nVertLevels)
        # We use advanced numpy indexing to grab the bottom layer for every cell
        num_cells = salinity_3d.shape[0]
        bottom_salinity = salinity_3d[np.arange(num_cells), max_level_cell]

        # Mask out land or areas where maxLevelCell was invalid (e.g., < 0)
        bottom_salinity = np.where(max_level_cell >= 0, bottom_salinity, np.nan)

        plot_data = bottom_salinity
        v_min = 34.0
        v_max = 35.0  # Typical Antarctic shelf salinity range, tweak as needed
        plotlabel = "Sea Floor Salinity [psu]"
        varcmap = "cmo.haline"  # cmocean haline is great for salinity
        norm = colors.Normalize(vmin=v_min, vmax=v_max)

    if opt_region == True:
        plot_data = np.where(iis == 0, np.nan, plot_data)
    # Save a clean version of bottom salinity for contouring later
    if var2plot == 'Sbot':
        salinity_contour_data = plot_data


    #PLOT
    projection = ccrs.PlateCarree()

    # define the transform that describes our dataset
    transform = ccrs.SouthPolarStereo()

    # create the figure and a GeoAxis
    fig, ax = plt.subplots(
        figsize=(10, 9),
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
        norm=norm,
        cmap=varcmap
    )

    if var2plot == 'Sbot':
        # Remove NaN values to prevent geometry errors in tricontour
        valid_mask = ~np.isnan(salinity_contour_data)

        if np.any(valid_mask):
            contour_line = ax.tricontour(
                lon[valid_mask],
                lat[valid_mask],
                salinity_contour_data[valid_mask],
                levels=[34.9],
                colors='black',
                linewidths=1.5,
                transform=ccrs.PlateCarree()  # Vital for mapping to your Cartopy axes properly
            )

            # Optional: Uncomment line below if you want inline text labels (e.g., "34.6")
            # ax.clabel(contour_line, inline=True, fontsize=8, fmt='%1.1f')

    # Because this is not a global mesh, it's necessary to explicitly set it's extent.
    ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())
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
    fig.colorbar(collection, fraction=0.1, shrink=0.5, label=plotlabel)
    if opt_save == 1:
        if opt_noGL == 1:
            plt.savefig(f'{fris_loc}/Fris_plots/snapshots_xy/{var2plot}_F{fnum}_Y{fy}_noGL.png', bbox_inches='tight',
                        dpi=300)
        else:
            plt.savefig(f'{fris_loc}/Fris_plots/snapshots_xy/{var2plot}_F{fnum}_Y{fy}.png', bbox_inches='tight',
                    dpi=300)
    else:
        plt.show()

