#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
#dataset = xr.open_dataset('/Users/irenavankova/Desktop/pyremap_test/temp_6000.0x5000.0km_10.0km_Antarctic_stereo_array.nc')
#dataset = xr.open_dataset('/Users/irenavankova/Desktop/pyremap_test/temp_oQU240.nc')
#dataset = xr.open_dataset('/Users/irenavankova/Desktop/pyremap_test/lifm_6000.0x5000.0km_10.0km_Antarctic_stereo.nc')

outdir = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_41-50_ts_1-50/custom2D_SOwISC12to60E2r4_to_6000.0x6000.0km_10.0km_Antarctic_stereo'
dataset = xr.open_dataset(f'{outdir}/mpaso_ANN_004101_005012_climo.nc')

lf = xr.open_dataset('/Users/irenavankova/Desktop/pyremap_test/test/lifm_6000.0x6000.0km_10.0km_Antarctic_stereo.nc')

landIceFloatingMask = numpy.squeeze(lf['landIceFloatingMask'].values)
#landIceFloatingMask[landIceFloatingMask < 1] = numpy.nan
x = dataset['x'].values
y = dataset['y'].values
#lat = dataset['lat'].values
#lon = dataset['lon'].values
#landIceFloatingMask = dataset['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature'][0, :, :].values  # Use the first time step (Time = 1)
Tb = dataset['timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature'].values  # Use the first time step (Time = 1)
ustar = dataset['timeMonthly_avg_landIceFrictionVelocity'].values  # Use the first time step (Time = 1)
Ti = dataset['timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature'].values  # Use the first time step (Time = 1)


cmap = plt.cm.seismic  # You can change this to another colormap if preferred
cmap.set_bad(color='gray')  # Set the 'bad' color (NaN) to gray

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(10, 8))


# Plot the landIceFloatingMask data
#c = ax.imshow((Tb-Ti)*landIceFloatingMask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
#              cmap=cmap, alpha=1, vmin=-0.5, vmax=0.5)
c = ax.imshow(landIceFloatingMask, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
              cmap=cmap, alpha=1, vmin=-0.5, vmax=0.5)

# Add a colorbar
fig.colorbar(c, ax=ax, label='U T*)')

# Title
#ax.set_title('Land Ice Floating Mask on Cartesian Coordinates')

# Display the plot
plt.show()