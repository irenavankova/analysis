#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

import matplotlib.colors as mcolors
import gmask_is

region_name = "Amundsen"
opt_save = 0
nsize = 0.1
nsize_sis = 0.5

# MPAS Ocean outputs
tclimo = ["006101_007012", "007101_008012"]
tseries = ["0001-0110", "0071-0080"]
tsegment = ["clim_61-70_ts_1-70", "clim_71-80_ts_1-80"]
sims = ["S12_control", "S12_control"]
#sims = ["S12_control"]


# Get region mask
if region_name == "Amery":
    sis_ctd = ["Amery_AM06"]
    iceshelves = ["Amery_shelf","Amery"]
    nshelf = 0
    y_lim = np.array([-2.5, 1.0])
    x_lim = np.array([33.0, 35.0])
elif region_name == "Ross":
    sis_ctd = []
    iceshelves = ["Ross_shelf","Ross"]
    nshelf = 0
    y_lim = np.array([-2.6, 1.6])
    x_lim = np.array([33.5, 35.2])
    nsize_sis = 0.1
elif region_name == "Filchner-Ronne":
    sis_ctd = ["Filchner_FNE1","Filchner_FNE3","Filchner_FSE1","Filchner_FSW1","Filchner_FSW2","Ronne_F1","Ronne_F2","Ronne_F3","Ronne_F4","Ronne_Site1","Ronne_Site2","Ronne_Site3","Ronne_Site4","Ronne_Site5"]
    iceshelves = ["Filchner-Ronne_shelf","Filchner-Ronne"]
    nshelf = 0
    y_lim = np.array([-2.7, 1.2])
    x_lim = np.array([33.5, 34.9])
    nsize_sis = 0.1
elif region_name == "Amundsen":
    sis_ctd = ["Thwaites_T1"]
    iceshelves = ["Amundsen_shelf","Amundsen"]
    nshelf = 0
    y_lim = np.array([-2.3, 2.0])
    x_lim = np.array([33.0, 35.2])
    nsize_sis = 0.5

iam, areaCell, isz = gmask_is.get_mask(iceshelves)

print(iam.shape)
print(isz.shape)

print(sum(np.sum(iam, axis=0)))
print(sum(np.sum(isz, axis=0)))

# MPAS pcean mesh
file_mesh = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'
dsMesh = xarray.open_dataset(file_mesh)
dsMesh = dsMesh[['nVertLevels','areaCell','maxLevelCell']]
dsMesh.load()
#areaCell = np.squeeze(dsMesh.areaCell.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)

# Sub ice shelf profile
for s in range(len(sis_ctd)):
    ctd_file = f'/Users/irenavankova/Work/data_sim/obs4E3SM/subshelf_ctd/nc_files/{sis_ctd[s]}.nc'
    dsCTD = xarray.open_dataset(ctd_file)
    dsCTD = dsCTD[['potentialTemperature','salinity']]
    dsCTD.load()
    PSsisN = np.squeeze(dsCTD.salinity.data)
    PTsisN = np.squeeze(dsCTD.potentialTemperature.data)
    if s == 0:
        PSsis = np.copy(PSsisN)
        PTsis = np.copy(PTsisN)
    else:
        PTsis = np.append(PTsis, PTsisN)
        PSsis = np.append(PSsis, PSsisN)

# Ship CTD data: Shenjie Zhou
sz_file = f'/Users/irenavankova/Desktop/Shenjie/Shenjie_CT.nc'
dsSZ = xarray.open_dataset(sz_file)
dsSZ = dsSZ[['ct']]
dsSZ.load()
CTsz = np.squeeze(dsSZ.ct.data)
sz_file = f'/Users/irenavankova/Desktop/Shenjie/Shenjie_SA.nc'
dsSZ = xarray.open_dataset(sz_file)
dsSZ = dsSZ[['sa']]
dsSZ.load()
SAsz = np.squeeze(dsSZ.sa.data)

SAsz = SAsz[:,np.squeeze(isz[nshelf,:])]
SAsz = SAsz.reshape((SAsz.shape[0]*SAsz.shape[1],)); SAsz = SAsz[~np.isnan(SAsz)]
CTsz = CTsz[:,np.squeeze(isz[nshelf,:])]
CTsz = CTsz.reshape((CTsz.shape[0]*CTsz.shape[1],)); CTsz = CTsz[~np.isnan(CTsz)]

PSsz = gsw.SP_from_SA(SAsz, p=0., lon=0., lat=-75.)
PTcz = gsw.pt_from_CT(SAsz, CTsz)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
plt.plot(PSsz, PTcz, color='gray', linestyle='None', marker='.', markersize=nsize)

# TS PLOT prep
PTbins = np.linspace(-2.8, 2, num=200)
PSbins = np.linspace(32.0, 35.5, num=200)
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

CS = plt.contour(PSgrid, PTgrid, neutralDensity, contours, linestyles=':', linewidths=0.5, colors='k', zorder=2)
plt.clabel(CS, fontsize=8, inline=1, fmt='%4.2f')

#clr = 'rybm'
clr = ["lightskyblue", "darkblue", "lightcoral", "maroon"]

ctr = 0
# Load and plot TS variables for control simulation(s)
for s in range(len(sims)):
    file_ts = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/mpaso_ANN_{tclimo[s]}_climo.nc'
    dsOutc = xarray.open_dataset(file_ts)
    dsOutc = dsOutc[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
    dsOutc.load()
    PTsim = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_temperature.data)
    PSsim = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_salinity.data)

    for n in range(len(iceshelves)):
        iis = iam[n,:]
        PTc = PTsim[iis]
        PSc = PSsim[iis]
        PTc = PTc.reshape((len(nVertLevels)*len(PTc),))
        PSc = PSc.reshape((len(nVertLevels)*len(PSc),))
        PTc = PTc[~np.isnan(PTc)]
        PSc = PSc[~np.isnan(PSc)]
        if n == nshelf:
            nsz = nsize
        else:
            nsz = nsize_sis
        plt.plot(PSc, PTc, clr[ctr], linestyle='None', marker='.', markersize=nsz)
        '''
        # Create the 2D histogram
        hist, xedges, yedges = np.histogram2d(PSc, PTc, bins=100)
        # Create the contour plot
        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
        # Plotting the contour lines
        plt.contour(X, Y, hist.T, levels=10, colors=clr[ctr])
        '''
        ctr = ctr + 1

#Plot subshelf CTDs
if 'PTsis' in globals():
    plt.plot(PSsis, PTsis, color='black', linestyle='None', marker='.', markersize=0.5)

plt.plot(PSbins, PTFreezing, linestyle='--', linewidth=1., color='g')
plt.ylim(y_lim)
plt.xlim(x_lim)
plt.show()

'''

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
    plt.xlim([33.8, 34.8])
elif reg == 'PIG':
    plt.ylim([-2, 1.2])
    plt.xlim([33.8, 34.8])
elif reg == 'Thwaites':
    plt.ylim([-2, 2])
    plt.xlim([33.8, 34.8])
elif reg == 'Ross':
    plt.ylim([-2.8, -1])
    plt.xlim([33.8, 35])
elif reg == 'FRIS':
    plt.ylim([-2.8, -1.5])
    plt.xlim([33.8, 35])
elif reg == 'Getz':
    plt.ylim([-2.5, 2])
    plt.xlim([33.8, 34.8])
else:
    plt.ylim([-2.5, 2])
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


'''


