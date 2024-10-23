#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

opt_save = 1
rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

t1 = 71
t1 = 61
t2 = t1+9
#trun = 'control'
#ttle = f'{trun} {t1}-{t2}'
#reg = 'Ross'
#reg = 'Amery'
#reg = 'PIG'
#reg = 'Thwaites'
reg = 'Totten'
ttle = f'{reg}_decade_{t1}_{t2}'
#p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/ECwISC30to60E2r1/ocean.ECwISC30to60E2r1.230220.nc'
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
#o_file = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_61-70_ts_1-70/mpaso_ANN_006101_007012_climo.nc'
o_file_c = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{t1}-{t2}_ts_1-{t2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'
o_file_m = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{t1}-{t2}_ts_1-{t2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'

#o_file = '/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_61-70_ts_1-70/mpaso_ANN_006101_007012_climo.nc'

dsMesh = xarray.open_dataset(p_file)
#dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','temperature','salinity','nVertLevels']]
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell']]
dsMesh.load()
areaCell = np.squeeze(dsMesh.areaCell.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)
maxLevelCell = np.squeeze(dsMesh.maxLevelCell.data)-1
#minLevelCell = np.squeeze(dsMesh.minLevelCell.data)
#print(maxLevelCell)

#Control
dsOutc = xarray.open_dataset(o_file_c)
dsOutc = dsOutc[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
dsOutc.load()
PTc = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_temperature.data)
PSc = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_salinity.data)
#Volume = np.squeeze(dsOut.timeMonthly_avg_layerThickness.data)
#Volume = Volume * np.transpose(np.array([areaCell,]*len(nVertLevels)))

#MALI
dsOutm = xarray.open_dataset(o_file_m)
dsOutm = dsOutm[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
dsOutm.load()
PTm = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_temperature.data)
PSm = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_salinity.data)
#Volume = np.squeeze(dsOut.timeMonthly_avg_layerThickness.data)
#Volume = Volume * np.transpose(np.array([areaCell,]*len(nVertLevels)))


#PT = np.squeeze(dsMesh.temperature.data)
#PS = np.squeeze(dsMesh.salinity.data)

lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
iii = (FloatingMask == 1)

if reg == 'Amery':
    iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
elif reg == 'PIG':
    iam = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880+360) & (lon < -98.4832+360)
elif reg == 'Thwaites':
    iam = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729+360) & (lon < -104.0573+360)
elif reg == 'Ross':
    iam = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212)
elif reg == 'FRIS':
    iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
elif reg == 'Totten':
    iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)

PTc = PTc[iam]
PSc = PSc[iam]
PTm = PTm[iam]
PSm = PSm[iam]
#Volume = Volume[iam]
maxLevelCell = maxLevelCell[iam]
#PT_floor = np.zeros((len(maxLevelCell)))
#PS_floor = np.zeros((len(maxLevelCell)))

#for v in range(len(maxLevelCell)):
#    PT_floor[v] = PT[v, maxLevelCell[v]]
#    PS_floor[v] = PS[v, maxLevelCell[v]]

# TS PLOT

PTc = PTc.reshape((len(nVertLevels)*len(PTc),))
PSc = PSc.reshape((len(nVertLevels)*len(PSc),))
PTm = PTm.reshape((len(nVertLevels)*len(PTm),))
PSm = PSm.reshape((len(nVertLevels)*len(PSm),))
#Volume = Volume.reshape((len(nVertLevels)*len(Volume),))

PTshelf = np.linspace(-2, np.nanmax(PTm), num=50)
PSshelf = np.linspace(-2, np.nanmax(PSm), num=50)
#iTmax = np.where(PTm==np.nanmax(PTm))
#PSshelf = np.linspace(33.8, PSm[iTmax], num=50)

PTbins = np.linspace(-2.8, 2, num=200)
PSbins = np.linspace(33.5, 35.2, num=200)

SAbins = gsw.SA_from_SP(PSbins, p=0., lon=0., lat=-75.)
CTbins = gsw.pt_from_CT(SAbins, PTbins)
CTgrid, SAgrid = numpy.meshgrid(CTbins, SAbins)

PSgrid = gsw.SP_from_SA(SAgrid, p=0., lon=0., lat=-75.)
PTgrid = gsw.pt_from_CT(SAgrid, CTgrid)

neutralDensity = gsw.sigma0(SAgrid, CTgrid)
rhoInterval = 0.2
contours = numpy.arange(23., 29.+rhoInterval, rhoInterval)

CTFreezing = gsw.CT_freezing(SAbins, 0, 1)
PTFreezing = gsw.t_from_CT(SAbins,CTFreezing, p=0.)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))

CS = plt.contour(PSgrid, PTgrid, neutralDensity, contours, linestyles=':', linewidths=0.5, colors='k', zorder=2)
plt.clabel(CS, fontsize=8, inline=1, fmt='%4.2f')

plt.plot(PSm, PTm, 'k.',markersize=0.5)
plt.plot(PSc, PTc, 'c.',markersize=0.5)
#plt.plot(PS_floor, PT_floor, 'r.')
#cmap = 'pink_r'
#hist, _, _, panel = plt.hist2d(PS, PT, bins=[PSbins, PTbins],weights=Volume, cmap=cmap, zorder=1)

plt.plot(PSbins, PTFreezing, linestyle='--', linewidth=1., color='m')
#plt.plot(PSshelf, PTshelf, linestyle='-', linewidth=1., color='c')

#plot GADE line (ambient melting)
Lw = 3.34 * pow(10,5)
cw = 3974.0
ci = 2009.0
TPi = -20

PTgade = PTshelf[-1]
PSgade = PSshelf[-1]
SAgade = gsw.SA_from_SP(PSgade, p=0., lon=0., lat=-75.)

CTgadeFreezing = gsw.CT_freezing(SAgade, 0, 1)
PTgadeFreezing = gsw.t_from_CT(SAgade,CTgadeFreezing, p=0.)

Teff = PTgade - PTgadeFreezing + Lw / cw + ci / cw * (PTgadeFreezing - TPi)
#Teff = 85
#Teff = 100;
print(Teff)

dTdS = 1 / PSgade * (Teff)
Tplot = PTgade + dTdS * (PSbins - PSgade)
plt.plot(PSbins, Tplot, linestyle='--', linewidth=1., color='b')

#plot SGR mixing line)

PTSGR = 0
PSSGR = 0
plt.plot([PSSGR, PSgade], [PTSGR, PTgade], linestyle='--', linewidth=1., color='r')

#plt.ylim([PTbins[0], PTbins[-1]])
#plt.xlim([PSbins[0], PSbins[-1]])

if reg == 'Amery':
    plt.ylim([-2.5, -0.8])
    plt.xlim([33.8, 34.5])
elif reg == 'PIG':
    plt.ylim([-2, 1.2])
    plt.xlim([33.8, 34.8])
elif reg == 'Thwaites':
    plt.ylim([-2, 2])
    plt.xlim([33.8, 34.8])
elif reg == 'Ross':
    plt.ylim([-2.8, 1])
    plt.xlim([33.8, 35])
elif reg == 'FRIS':
    plt.ylim([-2.8, -1.5])
    plt.xlim([33.8, 35])
else:
    plt.ylim([-2.5, 0])
    plt.xlim([33.8, 34.8])


fsize = 10
plt.xlabel('Salinity (PSU)', fontsize=fsize)
plt.ylabel('Potential temperature ($^\circ$C)', fontsize=fsize)
plt.title(ttle, fontsize=fsize)
plt.tight_layout()

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/global/TS_plots'
if opt_save == 1:

    plt.savefig(f'{dir_fig_save}/{ttle}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()





