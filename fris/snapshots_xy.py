#!/usr/bin/env python3

import xarray as xr
import numpy as np
import gmask_reg
import os

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# import cartopy.feature as cfeature

import cmocean
import mosaic
import matplotlib.path as mpath
from matplotlib import colors

fnum = 2
opt_noGL = 0
fy = '2-4'
opt_save = 1

# ADDED 'SpeedSurf' alongside 'SpeedBot' to variables to plot
variables_to_plot = ['SpeedSurf']

fris_loc = '/Users/ivankova/Desktop/Fris_hr'

mesh_file = f'/Users/ivankova/Desktop/Fris_hr/E3SM_init/ocean.ECwISC30to60E2r1.230220.nc'
out_file = mesh_file

opt_region = False

iceshelves = ["Shelf"]
if opt_region == True:
    iam = gmask_reg.get_mask(iceshelves, mesh_file, opt_noGL=opt_noGL, opt_wct=1)
    iis = iam[0, :]

dsMesh = xr.open_dataset(mesh_file)

dsMesh = dsMesh[
    ['latCell', 'lonCell', 'landIceFloatingMask', 'areaCell', 'maxLevelCell', 'nEdges', 'nVertices', 'lonEdge',
     'latEdge', 'lonVertex', 'latVertex', 'cellsOnEdge', 'cellsOnVertex', 'verticesOnEdge', 'verticesOnCell',
     'edgesOnVertex']]
dsMesh.load()

FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
areaCell = dsMesh['areaCell'].values
lat = dsMesh['latCell'].values
lon = dsMesh['lonCell'].values
lat = lat * 180 / np.pi
lon = lon * 180 / np.pi

max_level_cell = dsMesh['maxLevelCell'].values - 1  # Convert from 1-indexed to 0-indexed

ds = xr.open_dataset(out_file)

# Melt rate config
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
    elif var2plot == 'freeze':
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
        Tcb = ds['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature']
        Tci = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']
        Tc = (Tcb - Tci) * FloatingMask
        plot_data = Tc.isel(Time=0).values
        v_min = -1
        v_max = 1
        plotlabel = "Thermal driving [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Ti':
        Tci = ds['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature']
        plot_data = Tci.isel(Time=0).values
        v_min = -2.8
        v_max = -1.9
        plotlabel = "Temperature interface [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Tbl':
        Tcb = ds['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature']
        plot_data = Tcb.isel(Time=0).values
        v_min = -2.8
        v_max = -1.9
        plotlabel = "Temperature BL [C]"
        varcmap = "RdBu_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'ustar':
        Uc = ds['timeMonthly_avg_landIceFrictionVelocity']
        Uc = np.where(FloatingMask == 1, Uc * mtocm, np.nan)
        plot_data = Uc
        v_min = 0.4
        v_max = 1.2
        plotlabel = "Friction velocity [cm/s]"
        varcmap = "hot_r"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Sbot':
        if 'timeMonthly_avg_activeTracers_salinity' in ds:
            salinity_3d = ds['timeMonthly_avg_activeTracers_salinity'].isel(Time=0).values
        elif 'salinity' in ds:
            salinity_3d = ds['salinity'].isel(Time=0).values
        else:
            raise KeyError("Salinity variable not found.")

        num_cells = salinity_3d.shape[0]
        bottom_salinity = salinity_3d[np.arange(num_cells), max_level_cell]
        bottom_salinity = np.where(max_level_cell >= 0, bottom_salinity, np.nan)

        plot_data = bottom_salinity
        v_min = 34.0
        v_max = 35.0
        plotlabel = "Sea Floor Salinity [psu]"
        varcmap = "cmo.haline"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    elif var2plot == 'Tbot':
        if 'timeMonthly_avg_activeTracers_temperature' in ds:
            temp_3d = ds['timeMonthly_avg_activeTracers_temperature'].isel(Time=0).values
        elif 'temperature' in ds:
            temp_3d = ds['temperature'].isel(Time=0).values
        else:
            raise KeyError("Temperature variable not found.")

        num_cells = temp_3d.shape[0]
        bottom_temp = temp_3d[np.arange(num_cells), max_level_cell]
        bottom_temp = np.where(max_level_cell >= 0, bottom_temp, np.nan)

        plot_data = bottom_temp
        v_min = -2.6
        v_max = -1.0
        plotlabel = "Sea Floor Temp [C]"
        varcmap = "cmo.thermal"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)

    elif var2plot == 'SpeedBot':
        u_name = 'timeMonthly_avg_velocityZonal' if 'timeMonthly_avg_velocityZonal' in ds else 'velocityZonal'
        v_name = 'timeMonthly_avg_velocityMeridional' if 'timeMonthly_avg_velocityMeridional' in ds else 'velocityMeridional'

        if u_name not in ds or v_name not in ds:
            raise KeyError(f"Could not find velocity components '{u_name}' or '{v_name}' in the dataset.")

        u_3d = ds[u_name].isel(Time=0).values
        v_3d = ds[v_name].isel(Time=0).values

        num_cells = u_3d.shape[0]
        u_vector = u_3d[np.arange(num_cells), max_level_cell]
        v_vector = v_3d[np.arange(num_cells), max_level_cell]

        u_vector = np.where(max_level_cell >= 0, u_vector, np.nan)
        v_vector = np.where(max_level_cell >= 0, v_vector, np.nan)

        plot_data = np.sqrt(u_vector ** 2 + v_vector ** 2)
        v_min = 0.0
        v_max = 0.25
        plotlabel = "Bottom Flow Speed [m/s]"
        varcmap = "cmo.speed"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)

    # -------------------------------------------------------------------------
    # NEW CODE BLOCK: Surface Flow Speed Option
    # -------------------------------------------------------------------------
    elif var2plot == 'SpeedSurf':
        u_surf_name = 'surfaceVelocityZonal'
        v_surf_name = 'surfaceVelocityMeridional'

        if u_surf_name not in ds or v_surf_name not in ds:
            raise KeyError(
                f"Could not find surface velocity components '{u_surf_name}' or '{v_surf_name}' in the dataset.")

        # Since these are 2D (Time, nCells), we do not index by max_level_cell
        u_vector = ds[u_surf_name].isel(Time=0).values
        v_vector = ds[v_surf_name].isel(Time=0).values

        # Simple land masking check (using max_level_cell to drop zero-depth points)
        u_vector = np.where(max_level_cell >= 0, u_vector, np.nan)
        v_vector = np.where(max_level_cell >= 0, v_vector, np.nan)

        plot_data = np.sqrt(u_vector ** 2 + v_vector ** 2)
        v_min = 0.0
        v_max = 0.40  # Surface speeds are generally higher than bottom speeds; adjust dynamically
        plotlabel = "Surface Flow Speed [m/s]"
        varcmap = "cmo.speed"
        norm = colors.Normalize(vmin=v_min, vmax=v_max)
    # -------------------------------------------------------------------------

    if opt_region == True:
        plot_data = np.where(iis == 0, np.nan, plot_data)
    if var2plot == 'Sbot':
        salinity_contour_data = plot_data

    # PLOT Setup
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(
        figsize=(10, 9),
        constrained_layout=True,
        subplot_kw={"projection": projection},
    )

    descriptor = mosaic.Descriptor(dsMesh, projection=projection)

    collection = mosaic.polypcolor(
        ax,
        descriptor,
        plot_data,
        norm=norm,
        cmap=varcmap
    )

    if var2plot == 'Sbot':
        valid_mask = ~np.isnan(salinity_contour_data)
        if np.any(valid_mask):
            contour_line = ax.tricontour(
                lon[valid_mask],
                lat[valid_mask],
                salinity_contour_data[valid_mask],
                levels=[34.9],
                colors='black',
                linewidths=1.5,
                transform=ccrs.PlateCarree()
            )

        # -------------------------------------------------------------------------
        # FIXED UNIVERSAL VECTOR ARROW OVERLAP BLOCK
        # -------------------------------------------------------------------------
        if var2plot in ['SpeedBot', 'SpeedSurf']:
            # Filter within the Cartopy extent boundaries and drop NaNs
            spatial_mask = (lon >= -80) & (lon <= -25) & (lat >= -84) & (lat <= -70) & (~np.isnan(u_vector))

            # Stride density control (adjust based on mesh density)
            step = 50

            lon_q = lon[spatial_mask][::step]
            lat_q = lat[spatial_mask][::step]
            u_q = u_vector[spatial_mask][::step]
            v_q = v_vector[spatial_mask][::step]

            if len(lon_q) > 0:
                # CRITICAL FIX: Transform the lon/lat coordinates into map projection space first
                # This aligns the arrow anchor points to your exact figure view.
                xy_projected = projection.transform_points(ccrs.PlateCarree(), lon_q, lat_q)
                x_q = xy_projected[:, 0]
                y_q = xy_projected[:, 1]

                # Plot using the projected map coordinates instead of raw degrees
                q = ax.quiver(
                    x_q, y_q, u_q, v_q,
                    color='black',
                    scale=3.0,  # Decrease to 1.0 or 1.5 if arrows still look tiny
                    width=0.002,
                    headwidth=3,
                    pivot='middle'
                )
                # Add scale window reference
                ax.quiverkey(q, 0.85, 0.05, 0.2, r'0.2 m/s', labelpos='E', coordinates='axes', color='black')
                print(f"Successfully rendered {len(lon_q)} vector arrows for {var2plot}.")
            else:
                print(f"WARNING: No data points fell within the spatial mask window for {var2plot}.")
        # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    ax.set_extent([-80, -25, -84, -70], ccrs.PlateCarree())

    ax.set_aspect('auto')
    ax.gridlines(draw_labels=True)
    ax.set_facecolor('lightgray')

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