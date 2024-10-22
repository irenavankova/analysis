#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os


rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

#p_base = '/Users/irenavankova/Work/data_sim/SGR/test_fw_tracers/test_CB_FW_mpaso/'
#plot_folder = f'{p_base}/plots'
#fdir_sgr = f'{p_base}/sgr_on'
#fdir_ref = f'{p_base}/sgr_off'

opt_save = 1
Tbot = 1
dnum = '100'
#dir_now = f'f{dnum}/f{dnum}_00' #No melt no SGR
dir_now = f'f{dnum}/f{dnum}_010000B' #No melt yes SGR
#dir_now = f'f{dnum}/f{dnum}_10' #Yes melt no SGR
#dir_now = f'f{dnum}/f{dnum}_110000B' #Both


p_base = f'/Users/irenavankova/Work/data_sim/SGR/idealized/sg_data_conserve_04_yesC/{dir_now}'
dir_fig_save = f'/Users/irenavankova/Work/data_sim/SGR/idealized/plots/TS/{dir_now}/'
if not os.path.exists(dir_fig_save):
    os.makedirs(dir_fig_save)

ds = xarray.open_dataset(f'{p_base}/timeSeriesStatsMonthly.0002-12-01.nc')
ds.load()
dsMesh = xarray.open_dataset(f'{p_base}/restart.0003-01-01_00.00.00.nc')
dsMesh.load()
PT = np.squeeze(ds.timeMonthly_avg_activeTracers_temperature.data)
PS = np.squeeze(ds.timeMonthly_avg_activeTracers_salinity.data)
FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
nVertLevels = np.squeeze(dsMesh.nVertLevels.data)
print(nVertLevels)
#iii = np.where(FloatingMask == 15)
iii = (FloatingMask == 1)
PT = PT[iii,:]
print(PS.shape)
PS = PS[iii,:]
print(PS.shape)
print(len(PS))
PT = PT.reshape((len(nVertLevels)*len(PT),))
PS = PS.reshape((len(nVertLevels)*len(PS),))
print(PS.shape)


#PTshelf = np.linspace(-1.9, 1, num=50)
#PSshelf = np.linspace(33.8, 34.7, num=50)
PTshelf = np.linspace(-1.9, Tbot, num=50)
PSshelf = np.linspace(33.8, 34.7, num=50)

PTbins = np.linspace(-2, Tbot+0.2, num=50)
PSbins = np.linspace(33.6, 34.8, num=50)

SAbins = gsw.SA_from_SP(PSbins, p=0., lon=0., lat=-75.)
CTbins = gsw.pt_from_CT(SAbins, PTbins)
CTgrid, SAgrid = numpy.meshgrid(CTbins, SAbins)

#PSbins = gsw.SP_from_SA(SAbins, p=0., lon=0., lat=-75.)
#PSbins = gsw.SP_from_SA(SAbins, p=0., lon=0., lat=-75.)

PSgrid = gsw.SP_from_SA(SAgrid, p=0., lon=0., lat=-75.)
PTgrid = gsw.pt_from_CT(SAgrid, CTgrid)

#SA = gsw.SA_from_SP(PS, p=0., lon=0., lat=-75.)
#CT = gsw.CT_from_pt(SA, PT, p=0.)

neutralDensity = gsw.sigma0(SAgrid, CTgrid)
rhoInterval = 0.2
contours = numpy.arange(23., 29.+rhoInterval, rhoInterval)

CTFreezing = gsw.CT_freezing(SAbins, 0, 1)
PTFreezing = gsw.t_from_CT(SAbins,CTFreezing, p=0.)

fHeight = 5
fWidth = fHeight
plt.figure(figsize=(fWidth, fHeight))

CS = plt.contour(PSgrid, PTgrid, neutralDensity, contours, linewidths=1., colors='k', zorder=2)
plt.clabel(CS, fontsize=8, inline=1, fmt='%4.2f')

plt.plot(PS, PT, 'k.')

plt.plot(PSbins, PTFreezing, linestyle='--', linewidth=1., color='m')
plt.plot(PSshelf, PTshelf, linestyle='-', linewidth=1., color='c')

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



plt.ylim([PTbins[0], PTbins[-1]])
plt.xlim([PSbins[0], PSbins[-1]])

fsize = 8
plt.xlabel('Salinity (PSU)', fontsize=fsize)
plt.ylabel('Potential temperature ($^\circ$C)', fontsize=fsize)
plt.tight_layout()
if opt_save == 1:
    plt.savefig(f'{dir_fig_save}/TS.png', bbox_inches='tight', dpi=600)
else:
    plt.show()

'''

# do this before the inset because otherwise it moves the inset
# and cartopy doesn't play too well with tight_layout anyway

'''




