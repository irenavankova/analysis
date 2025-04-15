#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

import cartopy.crs as ccrs
import xarray as xr
import numpy
import cmocean

mtocm = 100
opt_save = 1
opt_mali = 1
opt_yy = 101

if opt_yy == 0:
    y1 = 21
    y2 = 21
    seas = '01'
    fname_in = f'mpaso_{seas}_00{y1}01_00{y2}01_climo.nc'
elif opt_yy == 22:
    y1 = 22
    y2 = 23
    seas = 'ANN'
    fname_in = f'mpaso_{seas}_00{y1}01_00{y2}12_climo.nc'
elif opt_yy == 41:
    y1 = 41
    y2 = 50
    seas = 'ANN'
    fname_in = f'mpaso_{seas}_00{y1}01_00{y2}12_climo.nc'
elif opt_yy == 101:
    y1 = 101
    y2 = 110
    seas = 'ANN'
    fname_in = f'mpaso_{seas}_0{y1}01_0{y2}12_climo.nc'

# Load the NetCDF file using xarray
lf = xr.open_dataset('/Users/irenavankova/Desktop/pyremap_test/test/lifm_6000.0x6000.0km_10.0km_Antarctic_stereo.nc')
landIceFloatingMask = numpy.squeeze(lf['landIceFloatingMask'].values)
print(sum(sum(landIceFloatingMask == 0)))
print(sum(sum(numpy.isnan(landIceFloatingMask))))
print(sum(sum(landIceFloatingMask == 1)))
landIceFloatingMask[(landIceFloatingMask > 0) & (landIceFloatingMask < 1)] = 1.0
binary_mask = numpy.copy(landIceFloatingMask)
landIceFloatingMask[landIceFloatingMask < 1] = numpy.nan

outdir_cont = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{y1}-{y2}_ts_1-{y2}/cc2D'
dc = xr.open_dataset(f'{outdir_cont}/{fname_in}')

if opt_mali == 8:
    mali_file = 'S12_mali_x8'
    fname = 'M8'
    cmax = 0.1
if opt_mali == 4:
    mali_file = 'S12_mali_x4'
    fname = 'M4'
    cmax = 0.08
if opt_mali == 1:
    mali_file = 'S12_mali'
    fname = 'M1'
    cmax = 0.02

outdir_mali = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{mali_file}/clim_{y1}-{y2}_ts_1-{y2}/cc2D'
dm = xr.open_dataset(f'{outdir_mali}/{fname_in}')

x = dc['x'].values
y = dc['y'].values

Tcb = dc['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature'].values  # Use the first time step (Time = 1)
Uc = dc['timeMonthly_avg_landIceFrictionVelocity'].values  # Use the first time step (Time = 1)
Uc = Uc*mtocm
Tci = dc['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature'].values  # Use the first time step (Time = 1)
Tc = Tcb-Tci

Tmb = dm['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature'].values  # Use the first time step (Time = 1)
Um = dm['timeMonthly_avg_landIceFrictionVelocity'].values  # Use the first time step (Time = 1)
Um = Um*mtocm
Tmi = dm['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature'].values  # Use the first time step (Time = 1)
Tm = Tmb-Tmi

Tp = Tm-Tc
Up = Um-Uc

#cmap = plt.cm.seismic

# Load the 'cmo.balance' colormap
cmapbal = cmocean.cm.balance
colors = cmapbal(numpy.linspace(0, 1, 256))
middle_idx = 128
colors[middle_idx] = numpy.array([1, 1, 1, 1])  # RGBA for white
cmap = plt.cm.colors.ListedColormap(colors)

#cmap = cmocean.cm.balance  # You can change this to another colormap if preferred
cmap = plt.cm.seismic


cmap.set_bad(color='green')  # Set the 'bad' color (NaN) to gray

# Set up the figure and axis
ncols = 2
nrows = 2
fHeight = 12
fWidth = 14
cm = 1/2.54

fig, ax = plt.subplots(nrows, ncols, figsize=(fWidth*cm, fHeight*cm))
plt.subplots_adjust(hspace=0.1*cm,wspace=0.1*cm)  # Increase vertical spacing
#gs = GridSpec(1, 1, figure=fig, hspace=-5*cm, wspace=-5*cm)

fsize = 8
ymin = 80
ymax = 530
xmin = 30
xmax = 580

# Define custom colormap with 3 colors: white for 0, gray for NaN, transparent for 1
colors = [(1, 1, 1),  # White for 0
          (0.5, 0.5, 0.5),  # Gray for NaN
          (1, 1, 1, 0)]  # Transparent for 1
# Create a new colormap from the defined colors
binary_cmap = ListedColormap(colors, name='binary_cmap')
binary_cmap.set_over((1, 1, 1, 0))  # Make values < 1 (NaN) transparent
binary_cmap.set_bad((0.5, 0.5, 0.5))  # Set NaN to gray

ctr = 0
for j in range(2):
    for k in range(2):
        if ctr == 0:
            ynow = Tc*Up
            ttl = '$u_p^{*}T_c^{*}$'
        elif ctr == 1:
            ynow = Tp*Uc
            ttl = '$u_c^{*}T_p^{*}$'
        elif ctr == 2:
            ynow = Tp*Up
            ttl = '$u_p^{*}T_p^{*}$'
        elif ctr == 3:
            ynow = (Tm * Um - Tc * Uc)
            ttl = '$u_s^{*}T_s^{*} - u_c^{*}T_c^{*}$'

        c = ax[j, k].imshow(ynow*landIceFloatingMask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
                            cmap=cmap, alpha=1, vmin=-cmax, vmax=cmax)
        ax[j,k].imshow(binary_mask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
                  cmap=binary_cmap, alpha=1)  # Adjust alpha for transparency
        ax[j, k].set_title(ttl, fontsize=fsize-1)
        ax[j, k].autoscale(enable=True, axis='both', tight=True)
        ax[j, k].tick_params(axis='both', labelsize=fsize)
        ax[j, k].set_xlim([xmin, xmax])
        ax[j, k].set_ylim([ymin, ymax])
        ctr = ctr + 1
        if k == 0:
            ax[j, k].set_ylabel('Northing (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_yticklabels(), visible=False)
        if j == nrows - 1:
            ax[j, k].set_xlabel('Easting (km)', fontsize=fsize)
        else:
            plt.setp(ax[j, k].get_xticklabels(), visible=False)


print(sum(sum(binary_mask == 0)))
print(sum(sum(numpy.isnan(binary_mask))))
print(sum(sum(binary_mask == 1)))

cbar = fig.colorbar(c, ax=ax, label='UT (cm/s * deg C)', fraction=0.02, pad=0.04)
cbar.ax.tick_params(labelsize=fsize-1)  # Set tick label font size
cbar.set_label('$u^{*}T^{*}$ ($^{\circ}$C m/s)', fontsize=fsize)  # Set colorbar label font size

if opt_save == 1:
    plt.savefig(f'/Users/irenavankova/Work/data_sim/SGR/global/Figs/melt_reynolds/melt_reynolds_{fname}_{y1}_{y2}_{seas}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()