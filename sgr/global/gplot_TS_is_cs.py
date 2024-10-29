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

opt_23 = 0

t1 = 71
t1 = 61
t2 = t1+9
ts2 = t2

if opt_23 == 1:
    t1 = 23
    t2 = 26
    ts2 = 32

r = 2

reglist = ["PIG", "Thwaites" , "Totten" , "Amery" , "Ross", "FRIS", "Getz", "Fimbul"]
reg = reglist[r]
ttle = f'{reg}_years_{t1}_{t2}'
p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
o_file_c = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_control/clim_{t1}-{t2}_ts_1-{ts2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'
o_file_m = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/S12_mali/clim_{t1}-{t2}_ts_1-{ts2}/mpaso_ANN_00{t1}01_00{t2}12_climo.nc'

dsMesh = xarray.open_dataset(p_file)
dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell','restingThickness']]
dsMesh.load()
areaCell = np.squeeze(dsMesh.areaCell.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)
H = np.squeeze(dsMesh.restingThickness.data)
H = np.nansum(H, axis=1)
H = np.squeeze(H)

#Masks
lat = np.squeeze(dsMesh.latCell.data)
lon = np.squeeze(dsMesh.lonCell.data)
lat = lat*180/np.pi
lon = lon*180/np.pi
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
iii = (FloatingMask == 1)

if reg == 'Amery':
    iam = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
    ics = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 70) & (lon < 80) & (H < 2000)
elif reg == 'PIG':
    iam = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880+360) & (lon < -98.4832+360)
    ics = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 252) & (lon < 256)  & (H < 2000)
elif reg == 'Thwaites':
    iam = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729+360) & (lon < -104.0573+360)
    ics = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 245) & (lon < 260)  & (H < 2000)
elif reg == 'Ross':
    #iam = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
    iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
    iam2 = (lat < -77.8) | (lon < 200)
    iam = np.logical_and(iam1, iam2)
    ics = (FloatingMask == 0) & (lat > -85.6) & (lat < -73) & (lon > 158.64) & (lon < 210) & (H < 2000)
elif reg == 'FRIS':
    iam = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
    ics = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (H < 2000)
elif reg == 'Totten':
    iam = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)
    ics = (FloatingMask == 0) & (lat > -68) & (lat < -65) & (lon > 114) & (lon < 124) & (H < 2000)
elif reg == 'Getz':
    iam = (FloatingMask == 1) & (lat > -75) & (lat < -73.5) & (lon > 225) & (lon < 245.5)
    #ics = (FloatingMask == 0) & (lat > -75) & (lat < -71) & (lon > 235) & (lon < 246) & (H < 2000)
    ics = (FloatingMask == 0) & (lat > -75.5) & (lat < -71) & (lon > 225) & (lon < 252)  & (H < 2000)
elif reg == 'Fimbul':
    iam = (FloatingMask == 1) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8))
    ics = (FloatingMask == 0) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8)) & (H < 2500)

# Start plot

PTbins = np.linspace(-2.8, 2.5, num=200)
PSbins = np.linspace(33.2, 35.2, num=200)

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

# Add field values
for s in range(2):
    for v in range(2):
        if v == 1:
            # Control
            dsOutc = xarray.open_dataset(o_file_c)
            dsOutc = dsOutc[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
            dsOutc.load()
            PT = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_temperature.data)
            PS = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_salinity.data)
            clr = 'm'
        else:
            #MALI
            dsOutm = xarray.open_dataset(o_file_m)
            dsOutm = dsOutm[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
            dsOutm.load()
            PT = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_temperature.data)
            PS = np.squeeze(dsOutm.timeMonthly_avg_activeTracers_salinity.data)
            clr = 'k'


        if s == 1:
            inow = iam
            if clr == 'm':
                clr = 'c'
            else:
                clr = 'b'
        else:
            inow = ics
        PTc = PT[inow]
        PSc = PS[inow]

# TS PLOT

        PTc = PTc.reshape((len(nVertLevels)*len(PTc),))
        PSc = PSc.reshape((len(nVertLevels)*len(PSc),))

        if (v == 0) & (s == 1):
            PTshelf = np.linspace(-2, np.nanmax(PTc), num=50)
            PSshelf = np.linspace(-2, np.nanmax(PSc), num=50)

        plt.plot(PSc, PTc, f'{clr}', linestyle='', marker='.', markersize=0.5)

#plot GADE line (ambient melting)
Lw = 3.34 * pow(10,5)
cw = 3974.0
ci = 2009.0
TPi = -20

print(PTshelf[-1])
PTgade = PTshelf[-1]
PSgade = PSshelf[-1]
SAgade = gsw.SA_from_SP(PSgade, p=0., lon=0., lat=-75.)
CTgadeFreezing = gsw.CT_freezing(SAgade, 0, 1)
PTgadeFreezing = gsw.t_from_CT(SAgade,CTgadeFreezing, p=0.)
Teff = PTgade - PTgadeFreezing + Lw / cw + ci / cw * (PTgadeFreezing - TPi)
#Teff = 85
#Teff = 100;
dTdS = 1 / PSgade * (Teff)
Tplot = PTgade + dTdS * (PSbins - PSgade)
plt.plot(PSbins, Tplot, linestyle=':', linewidth=1., color='g')
#plot SGR mixing line)

PTSGR = 0
PSSGR = 0
plt.plot([PSSGR, PSgade], [PTSGR, PTgade], linestyle=':', linewidth=1., color='g')

plt.plot(PSbins, PTFreezing, linestyle=':', linewidth=1., color='g')

if reg == 'Amery':
    plt.ylim([-2.5, 1.5])
    plt.xlim([33.5, 34.8])
elif reg == 'PIG':
    plt.ylim([0, 2])
    plt.xlim([34.4, 34.8])
elif reg == 'Thwaites':
    plt.ylim([-2, 2])
    plt.xlim([33.8, 34.8])
elif reg == 'Ross':
    plt.ylim([-2.8, 1.8])
    plt.xlim([33.6, 35])
elif reg == 'FRIS':
    plt.ylim([-2.8, 1.5])
    plt.xlim([33.6, 35])
elif reg == 'Getz':
    plt.ylim([0.5, 2])
    plt.xlim([34.4, 34.8])
elif reg == 'Fimbul':
    plt.ylim([-2.5, 2])
    plt.xlim([33.5, 34.8])
elif reg == 'Totten':
    plt.ylim([-2.5, 1])
    plt.xlim([33.7, 34.8])
else:
    plt.ylim([-2.5, 2])
    plt.xlim([33.5, 34.8])

fsize = 10
plt.xlabel('Salinity (PSU)', fontsize=fsize)
plt.ylabel('Potential temperature ($^\circ$C)', fontsize=fsize)
plt.title(ttle, fontsize=fsize)
plt.tight_layout()

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/global/TS_plots_cs'
if opt_save == 1:

    plt.savefig(f'{dir_fig_save}/TS_cs_{ttle}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()





