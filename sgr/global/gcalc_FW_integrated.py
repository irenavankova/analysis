#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import cmocean

chic = 1
opt_mali = 1
N = 13

if opt_mali == 8:
    mali_file = 'S12_mali_x8'
    rname = '20241029.GMPAS-JRA1p5-DIB-PISMF-DSGR.TL319_SOwISC12to60E2r4.sgr_mali_x8.chicoma-cpu'
    fname = 'M8'
    N = 30 * 12
if opt_mali == 4:
    mali_file = 'S12_mali_x4'
    rname = '20241022.GMPAS-JRA1p5-DIB-PISMF-DSGR.TL319_SOwISC12to60E2r4.sgr_mali_x4.chicoma-cpu'
    fname = 'M4'
    N = 30 * 12
if opt_mali == 1:
    mali_file = 'S12_mali'
    rname = '20240802.GMPAS-JRA1p5-DIB-PISMF-DSGR.TL319_SOwISC12to60E2r4.sgr_mali.chicoma-cpu'
    fname = 'M1'
    N = 90 * 12
if opt_mali == 0:
    mali_file = 'S12_control'
    rname = '20240802.GMPAS-JRA1p5-DIB-PISMF.TL319_SOwISC12to60E2r4.sgr_control.chicoma-cpu'
    fname = 'M0'
    N = 90 * 12

# MPAS Ocean mesh file
if chic == 1:
    cdir = f'/lustre/scratch5/vankova/E3SM/scratch/chicoma-cpu/{rname}'
    p_file = f'{cdir}/run/{rname}.mpaso.rst.0021-01-01_00000.nc'
    cpref = f'{rname}.mpaso.hist.am.timeSeriesStatsMonthly.'
    sdir = '/lustre/scratch5/vankova/E3SM/scratch/plots/sgr_fw'
    csave = f'{sdir}/fw_tseries_{fname}.nc'
else:
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
    csave = f'/Users/irenavankova/Desktop/bla/bla.nc'
    N = 2

time = np.arange(N)

dsMesh = xr.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','areaCell']]
dsMesh.load()

lat = np.squeeze(dsMesh['latCell'].values)
lon = np.squeeze(dsMesh['lonCell'].values)
areaCell = np.squeeze(dsMesh['areaCell'].values)
FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
lat = lat*180/np.pi
lon = lon*180/np.pi

VFWLIicav = np.zeros([N])
VFWLIishelf = np.zeros([N])
VFWLIis60 = np.zeros([N])
VFWLI500is60 = np.zeros([N])
VFWLI1000is60 = np.zeros([N])
VFWSRicav = np.zeros([N])
VFWSRishelf = np.zeros([N])
VFWSRis60 = np.zeros([N])
VFWSR500is60 = np.zeros([N])
VFWSR1000is60 = np.zeros([N])

# FW Ant projected file
yy = 21
mm = 1
for j in range(N):
    print(j)
    if chic == 1:
        p_file = f'{cdir}/run/{cpref}{str(yy).zfill(4)}-{str(mm).zfill(2)}-01.nc'
        mm = mm + 1
        if mm == 13:
            yy = yy + 1
            mm = 1

    else:
        if j == 0:
            y1 = '41'
            y2 = '50'
        else:
            y1 = '101'
            y2 = '110'
        fdir = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{y1}-{y2}_ts_1-{y2}/fw'
        p_file = f'{fdir}/mpaso_ANN_{str(y1).zfill(4)}01_{str(y2).zfill(4)}12_climo.nc'

    print(p_file)
    ds = xr.open_dataset(p_file)
    layerThicknessC = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)

    if j == 0:
        sshC = np.squeeze(ds['timeMonthly_avg_ssh'].values)
        depthC = layerThicknessC
        depthC[:,0] = depthC[:,0] - sshC
        depthC = np.cumsum(depthC, axis=1)
        depthC = np.nanmax(depthC, axis=1)

        icav = (FloatingMask == 1)
        ishelf = ((depthC < 1500) & (lat < -60.0)) | (FloatingMask == 1)
        is60 = (lat < -60.0)

        '''
        fHeight = 5
        fWidth = fHeight
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[icav],lat[icav],'r.')

        #South of 60 and cavity or shallower than 1500m
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[ishelf],lat[ishelf],'g.')

        #South of 60 and cavity or shallower than 1500m
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[is60],lat[is60],'g.')
        #plt.show()
        '''

    #load
    lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
    srfwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_subglacialRunoffFreshWaterConcentration'].values)
    layerThicknessC = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)

    #adjust bad values
    lifwc[lifwc < 0] = 0.0
    lifwc = np.nan_to_num(lifwc, nan=0.0)
    srfwc[srfwc < 0] = 0.0
    srfwc = np.nan_to_num(srfwc, nan=0.0)

    #Calc LFW thickness over different regions
    HFWLI = np.nansum(np.multiply(layerThicknessC, lifwc), axis=1)
    HFWLI500 = np.nansum(np.multiply(layerThicknessC[:,0:32], lifwc[:,0:32]), axis=1)
    HFWLI1000 = np.nansum(np.multiply(layerThicknessC[:,0:39], lifwc[:,0:39]), axis=1)

    VFWLIicav[j] = np.nansum(np.multiply(HFWLI[icav],areaCell[icav]))
    VFWLIishelf[j] = np.nansum(np.multiply(HFWLI[ishelf],areaCell[ishelf]))
    VFWLIis60[j] = np.nansum(np.multiply(HFWLI[is60],areaCell[is60]))
    VFWLI500is60[j] = np.nansum(np.multiply(HFWLI500[is60],areaCell[is60]))
    VFWLI1000is60[j] = np.nansum(np.multiply(HFWLI1000[is60],areaCell[is60]))

    # Calc SFW thickness over different regions
    HFWSR = np.nansum(np.multiply(layerThicknessC, srfwc), axis=1)
    HFWSR500 = np.nansum(np.multiply(layerThicknessC[:, 0:32], srfwc[:, 0:32]), axis=1)
    HFWSR1000 = np.nansum(np.multiply(layerThicknessC[:, 0:39], srfwc[:, 0:39]), axis=1)

    VFWSRicav[j] = np.nansum(np.multiply(HFWSR[icav], areaCell[icav]))
    VFWSRishelf[j] = np.nansum(np.multiply(HFWSR[ishelf], areaCell[ishelf]))
    VFWSRis60[j] = np.nansum(np.multiply(HFWSR[is60], areaCell[is60]))
    VFWSR500is60[j] = np.nansum(np.multiply(HFWSR500[is60], areaCell[is60]))
    VFWSR1000is60[j] = np.nansum(np.multiply(HFWSR1000[is60], areaCell[is60]))

plt.plot(time,VFWLIicav,label = 'VFWLIicav')
plt.plot(time,VFWLIishelf,label = 'VFWLIishelf')
plt.plot(time,VFWLIis60,label = 'VFWLIis60')
plt.plot(time,VFWLI500is60,label = 'VFWLI500is60')
plt.plot(time,VFWLI1000is60,label = 'VFWLI1000is60')

plt.plot(time,VFWSRicav,'--',label = 'VFWSRicav')
plt.plot(time,VFWSRishelf,'--',label = 'VFWSRishelf')
plt.plot(time,VFWSRis60,'--',label = 'VFWSRis60')
plt.plot(time,VFWSR500is60,'--',label = 'VFWSR500is60')
plt.plot(time,VFWSR1000is60,'--',label = 'VFWSR1000is60')

plt.legend(loc=2, prop={'size': 8})
if chic == 1:
    plt.savefig(f'{sdir}/fw_tseries_plot_{fname}.png',bbox_inches='tight', dpi=600)
else:
    plt.show()

print('save to netcdf')
d1 = xr.DataArray(VFWLIicav, dims=("time"),
                                coords={'time': time},
                                name='VFWLIicav')  # Set the variable name
d2 = xr.DataArray(VFWLIishelf, dims=("time"),
                                coords={'time': time},
                                name='VFWLIishelf')  # Set the variable name
d3 = xr.DataArray(VFWLIis60, dims=("time"),
                                coords={'time': time},
                                name='VFWLIis60')  # Set the variable name
d4 = xr.DataArray(VFWLI500is60, dims=("time"),
                                coords={'time': time},
                                name='VFWLI500is60')  # Set the variable name
d5 = xr.DataArray(VFWLI1000is60, dims=("time"),
                                coords={'time': time},
                                name='VFWLI1000is60')  # Set the variable name
d6 = xr.DataArray(VFWSRicav, dims=("time"),
                                coords={'time': time},
                                name='VFWSRicav')  # Set the variable name
d7 = xr.DataArray(VFWSRishelf, dims=("time"),
                                coords={'time': time},
                                name='VFWSRishelf')  # Set the variable name
d8 = xr.DataArray(VFWSRis60, dims=("time"),
                                coords={'time': time},
                                name='VFWSRis60')  # Set the variable name
d9 = xr.DataArray(VFWSR500is60, dims=("time"),
                                coords={'time': time},
                                name='VFWSR500is60')  # Set the variable name
d10 = xr.DataArray(VFWSR1000is60, dims=("time"),
                                coords={'time': time},
                                name='VFWSR1000is60')  # Set the variable name

ds = xr.Dataset({
    'VFWLIicav': d1,
    'VFWLIishelf': d2,
    'VFWLIis60': d3,
    'VFWLI500is60': d4,
    'VFWLI1000is60': d5,
    'VFWSRicav': d6,
    'VFWSRishelf': d7,
    'VFWSRis60': d8,
    'VFWSR500is60': d9,
    'VFWSR1000is60': d10
})

ds.to_netcdf(f"{csave}")

#print('Looping')
#for i in range(FWC.shape[0]):