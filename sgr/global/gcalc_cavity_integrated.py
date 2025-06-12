#!/usr/bin/env python3

# Load the NetCDF file using xarray
# Replace 'your_file.nc' with the actual path to your NetCDF file
import matplotlib.pyplot as plt

import xarray as xr
import numpy as np

chic = 1
opt_mali = 8
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
    sdir = '/lustre/scratch5/vankova/E3SM/scratch/plots/sgr_si'
    csave = f'{sdir}/si_cav_tseries_{fname}.nc'
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

VFWLI_Amery = np.zeros([N])
VFWLI_Ross = np.zeros([N])
VFWLI_FRIS = np.zeros([N])
VFWSI_Amery = np.zeros([N])
VFWSI_Ross = np.zeros([N])
VFWSI_FRIS = np.zeros([N])

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
            y1 = '50'
            y2 = '50'
        else:
            y1 = '110'
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

        iAmery = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
        iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
        iam2 = (lat < -77.8) | (lon < 200)
        iRoss = np.logical_and(iam1, iam2)
        iFRIS = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)

        '''
        fHeight = 5
        fWidth = fHeight
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[iAmery],lat[iAmery],'r.')

        #South of 60 and cavity or shallower than 1500m
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[iRoss],lat[iRoss],'g.')

        #South of 60 and cavity or shallower than 1500m
        plt.figure(figsize=(fWidth, fHeight))
        plt.plot(lon,lat,'b.')
        plt.plot(lon[iFRIS],lat[iFRIS],'g.')
        #plt.show()
        '''


    #load
    lifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration'].values)
    sifwc = np.squeeze(ds['timeMonthly_avg_freshwaterTracers_seaIceFreshWaterConcentration'].values)
    layerThicknessC = np.squeeze(ds['timeMonthly_avg_layerThickness'].values)

    #adjust bad values
    lifwc[lifwc < 0] = 0.0
    lifwc = np.nan_to_num(lifwc, nan=0.0)
    sifwc[sifwc < -1.0] = 0.0
    sifwc[sifwc > 1.0] = 0.0
    sifwc = np.nan_to_num(sifwc, nan=0.0)

    #Calc LFW thickness over different regions
    HFWLI = np.nansum(np.multiply(layerThicknessC, lifwc), axis=1)
    HFWSI = np.nansum(np.multiply(layerThicknessC, sifwc), axis=1)

    VFWLI_Amery[j] = np.nansum(np.multiply(HFWLI[iAmery],areaCell[iAmery]))
    VFWLI_Ross[j] = np.nansum(np.multiply(HFWLI[iRoss],areaCell[iRoss]))
    VFWLI_FRIS[j] = np.nansum(np.multiply(HFWLI[iFRIS],areaCell[iFRIS]))

    VFWSI_Amery[j] = np.nansum(np.multiply(HFWSI[iAmery], areaCell[iAmery]))
    VFWSI_Ross[j] = np.nansum(np.multiply(HFWSI[iRoss], areaCell[iRoss]))
    VFWSI_FRIS[j] = np.nansum(np.multiply(HFWSI[iFRIS], areaCell[iFRIS]))

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(time,VFWLI_Amery,label = 'VFWLI_Amery')
plt.plot(time,VFWLI_Ross,label = 'VFWLI_Ross')
plt.plot(time,VFWLI_FRIS,label = 'VFWLI_FRIS')

plt.plot(time,VFWSI_Amery,label = 'VFWSI_Amery')
plt.plot(time,VFWSI_Ross,label = 'VFWSI_Ross')
plt.plot(time,VFWSI_FRIS,label = 'VFWSI_FRIS')

plt.legend(loc=2, prop={'size': 8})
if chic == 1:
    plt.savefig(f'{sdir}/si_cav_tseries_plot_{fname}.png',bbox_inches='tight', dpi=600)
else:
    plt.show()

print('save to netcdf')
d1 = xr.DataArray(VFWLI_Amery, dims=("time"),
                                coords={'time': time},
                                name='VFWLI_Amery')  # Set the variable name
d2 = xr.DataArray(VFWLI_Ross, dims=("time"),
                                coords={'time': time},
                                name='VFWLI_Ross')  # Set the variable name
d3 = xr.DataArray(VFWLI_FRIS, dims=("time"),
                                coords={'time': time},
                                name='VFWLI_FRIS')  # Set the variable name
d4 = xr.DataArray(VFWSI_Amery, dims=("time"),
                                coords={'time': time},
                                name='VFWSI_Amery')  # Set the variable name
d5 = xr.DataArray(VFWSI_Ross, dims=("time"),
                                coords={'time': time},
                                name='VFWSI_Ross')  # Set the variable name
d6 = xr.DataArray(VFWSI_FRIS, dims=("time"),
                                coords={'time': time},
                                name='VFWSI_FRIS')  # Set the variable name


ds = xr.Dataset({
    'VFWLI_Amery': d1,
    'VFWLI_Ross': d2,
    'VFWLI_FRIS': d3,
    'VFWSI_Amery': d4,
    'VFWSI_Ross': d5,
    'VFWSI_FRIS': d6
})

ds.to_netcdf(f"{csave}")

#print('Looping')
#for i in range(FWC.shape[0]):