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

region_name = "Fimbul"
appe = "_sgr_41_50"
opt_save = 1
nsize = 0.1
nsize_sis = 0.5

# MPAS Ocean outputs

c1 = np.array([41, 41, 41, 41])
#c1 = np.array([21, 41, 101])
c2 = c1 + 9
t2 = c1 + 9

tseries = []
tsegment = []
tclimo = []
for n in range(len(c1)):
    tseries.append(f"0001-{t2[n]:04}")
    tsegment.append(f"clim_{c1[n]}-{c2[n]}_ts_1-{t2[n]}")
    tclimo.append(f"{c1[n]:04}01_{c2[n]:04}12")

#tclimo = ["006101_007012", "007101_008012"]
#tseries = ["0001-0110", "0071-0080"]
#tsegment = ["clim_61-70_ts_1-70", "clim_71-80_ts_1-80"]
sims = ["S12_control", "S12_mali", "S12_mali_x4", "S12_mali_x8"]
#sims = ["S12_control"]


# Get region mask
if region_name == "Amery":
    ttl = "Amery"
    sis_ctd = ["Amery_AM06"]
    iceshelves = ["Amery_shelf","Amery"]
    nshelf = 0
    y_lim = np.array([-2.5, 1.0])
    x_lim = np.array([33.75, 34.75])
    nsize = 0.5
elif region_name == "Ross":
    ttl = "Ross"
    sis_ctd = []
    iceshelves = ["Ross_shelf","Ross"]
    nshelf = 0
    y_lim = np.array([-2.6, 1.6])
    x_lim = np.array([33.8, 35.0])
    nsize = 0.2
    nsize_sis = 0.2
elif region_name == "FRIS":
    ttl = "Filchner-Ronne"
    sis_ctd = ["Filchner_FNE1","Filchner_FNE3","Filchner_FSE1","Filchner_FSW1","Filchner_FSW2","Ronne_F1","Ronne_F2","Ronne_F3","Ronne_F4","Ronne_Site1","Ronne_Site2","Ronne_Site3","Ronne_Site4","Ronne_Site5"]
    iceshelves = ["Filchner-Ronne_shelf","Filchner-Ronne"]
    nshelf = 0
    y_lim = np.array([-2.7, 1.2])
    x_lim = np.array([33.6, 34.9])
    nsize = 0.2
    nsize_sis = 0.2
elif region_name == "Amundsen":
    ttl = "Amundsen"
    sis_ctd = ["Thwaites_T1"]
    iceshelves = ["Amundsen_shelf","Amundsen"]
    nshelf = 0
    y_lim = np.array([-2.3, 2.0])
    x_lim = np.array([33.5, 34.75])
    nsize_sis = 0.5
elif region_name == "Larsen":
    ttl = "Larsen C"
    sis_ctd = ["LarsenC_F111","LarsenC_F211"]
    iceshelves = ["Larsen_C_shelf","Larsen_C"]
    nshelf = 0
    y_lim = np.array([-2.2, 0.5])
    x_lim = np.array([33.7, 34.8])
    nsize_sis = 0.1
    nsize = 0.3
elif region_name == "TottenMU":
    ttl = "Totten/Moscow University"
    sis_ctd = []
    iceshelves = ["TottenMU_shelf","TottenMU"]
    nshelf = 0
    y_lim = np.array([-2.1, 0.6])
    x_lim = np.array([33.7, 34.8])
    nsize_sis = 0.5
    nsize = 0.3
elif region_name == "Fimbul":
    ttl = "Fimbul"
    sis_ctd = ["Fimbul_M1","Fimbul_M2","Fimbul_M3"]
    iceshelves = ["Fimbul_shelf","Fimbul"]
    nshelf = 0
    y_lim = np.array([-2.1, 1.3])
    x_lim = np.array([33.7, 34.8])
    nsize_sis = 0.3
    nsize = 0.3

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
areaCell = np.squeeze(dsMesh.areaCell.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)

# Sub ice shelf profile
'''
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
'''

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))
if 'PSsz' in globals():
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
#clr = ["lightskyblue", "royalblue", "moccasin", "darkorange", "yellowgreen","darkolivegreen","plum", "purple", "lightcoral", "maroon"]
clr = ["lightcoral", "brown", "moccasin", "darkorange", "lightskyblue", "dodgerblue", "plum", "indigo"]

ctr = 0
PTcm = np.zeros(len(tseries)*len(iceshelves))
PScm = np.zeros(len(tseries)*len(iceshelves))
PTcmn = np.zeros(len(tseries))
PScmn = np.zeros(len(tseries))
PTcms = np.zeros(len(tseries))
PScms = np.zeros(len(tseries))
# Load and plot TS variables for control simulation(s)
for s in range(len(tseries)):
    file_ts = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/SGR/ncfiles/{sims[s]}/{tsegment[s]}/mpaso_ANN_{tclimo[s]}_climo.nc'
    dsOutc = xarray.open_dataset(file_ts)
    dsOutc = dsOutc[['timeMonthly_avg_activeTracers_temperature','timeMonthly_avg_activeTracers_salinity','timeMonthly_avg_layerThickness']]
    dsOutc.load()
    PTsim = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_temperature.data)
    PSsim = np.squeeze(dsOutc.timeMonthly_avg_activeTracers_salinity.data)
    Volume = np.squeeze(dsOutc.timeMonthly_avg_layerThickness.data)
    Volume = Volume * np.transpose(np.array([areaCell,]*len(nVertLevels)))

    for n in range(len(iceshelves)):
        iis = iam[n,:]
        PTc = PTsim[iis]
        PSc = PSsim[iis]
        Volc = Volume[iis]
        PTc = PTc.reshape((len(nVertLevels)*len(PTc),))
        PSc = PSc.reshape((len(nVertLevels)*len(PSc),))
        Volc = Volc.reshape((len(nVertLevels) * len(Volc),))
        PTc = PTc[~np.isnan(PTc)]
        PSc = PSc[~np.isnan(PSc)]
        Volc = Volc[~np.isnan(Volc)]
        PTcm[ctr] = np.dot(PTc, Volc) / np.sum(Volc)
        PScm[ctr] = np.dot(PSc, Volc) / np.sum(Volc)
        if n == nshelf:
            nsz = nsize
            PTcmn[s] = PTcm[ctr]
            PScmn[s] = PScm[ctr]
        else:
            nsz = nsize_sis
            PTcms[s] = PTcm[ctr]
            PScms[s] = PScm[ctr]
        plt.plot(PSc, PTc, clr[ctr], linestyle='None', marker='.', markersize=nsz)
        #plt.plot(PScm, PTcm, clr[ctr], linestyle='None', marker='s', markersize=8, mfc=clr[ctr], mec='k')
        '''
        # Create the 2D histogram
        hist, xedges, yedges = np.histogram2d(PSc, PTc, bins=100)
        # Create the contour plot
        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
        # Plotting the contour lines
        plt.contour(X, Y, hist.T, levels=10, colors=clr[ctr])
        '''
        ctr = ctr + 1

ctr = 0
plt.plot(PScms, PTcms, 'k')
plt.plot(PScmn, PTcmn, 'k')
for s in range(len(tseries)):
    for n in range(len(iceshelves)):
        if n == nshelf:
            mrk = 's'
        else:
            mrk = 'o'
        plt.plot(PScm[ctr], PTcm[ctr], clr[ctr], linestyle='None', marker=mrk, markersize=5, mfc=clr[ctr], mec='k')
        ctr = ctr + 1
#Plot subshelf CTDs
if 'PTsis' in globals():
    plt.plot(PSsis, PTsis, color='black', linestyle='None', marker='.', markersize=0.5)

plt.plot(PSbins, PTFreezing, linestyle='--', linewidth=1., color='g')
plt.ylim(y_lim)
plt.xlim(x_lim)

fsize = 10
plt.xlabel('Salinity (PSU)', fontsize=fsize+4)
plt.ylabel('Potential temperature ($^\circ$C)', fontsize=fsize+4)
plt.title(ttl, fontsize=fsize+4)
plt.tight_layout()

dir_fig_save = '/Users/irenavankova/Work/data_sim/SGR/global/TS_eval_sgr'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/{region_name}{appe}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

