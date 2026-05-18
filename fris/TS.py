#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

#import sys
#sys.path.append('../sgr/global')
import gmask_reg

region_name = "FRIS"
opt_save = 1
nsize = 0.1
nsize_sis = 0.5
nsize_obs = 0.1

# MPAS Ocean outputs

tsegment = ["clim_2005-2014_ts_1951-2014", "clim_2091-2100_ts_2015-2100", "clim_2091-2100_ts_2015-2100"]
sims = ["hist", "fismf_701", "pismf"]
appe = "_ssp_91_100"
#sims = ["S12_control"]

# Get region mask
if region_name == "FRIS":
    ttl = "Filchner-Ronne"
    sis_ctd = ["Filchner_FNE1","Filchner_FNE3","Filchner_FSE1","Filchner_FSW1","Filchner_FSW2","Ronne_F1","Ronne_F2","Ronne_F3","Ronne_F4","Ronne_Site1","Ronne_Site2","Ronne_Site3","Ronne_Site4","Ronne_Site5"]
    iceshelves = ["Filchner-Ronne_shelf","Filchner-Ronne"]
    nshelf = 0
    y_lim = np.array([-2.7, 2.2])
    x_lim = np.array([33.5, 34.9])
    nsize = 0.2
    nsize_sis = 0.2


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

# Plot info

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))

# TS PLOT prep
PTbins = np.linspace(-2.8, 4, num=200)
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
PTcm = np.zeros(len(sims)*len(iceshelves))
PScm = np.zeros(len(sims)*len(iceshelves))
PTcmn = np.zeros(len(sims))
PScmn = np.zeros(len(sims))
PTcms = np.zeros(len(sims))
PScms = np.zeros(len(sims))
# Load and plot TS variables for control simulation(s)
for s in range(len(sims)):
    file_ts = f'/Users/irenavankova/Work/data_sim/E3SM_outputs/FISMF/ncfiles/{sims[s]}/{tsegment[s]}/TS.nc'
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
#plt.plot(PScms, PTcms, 'k')
#plt.plot(PScmn, PTcmn, 'k')
for s in range(len(sims)):
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

dir_fig_save = '/Users/irenavankova/Work/data_sim/FISMF/TS_sim'
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/{region_name}{appe}.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

