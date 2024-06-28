#!/usr/bin/env python3
import numpy as np
import xarray
import os
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt

rho_fw = 1000.
secPerYear = 365 * 24 * 60 * 60

p_base = '/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC'
sgr_base = '/Users/irenavankova/Work/data_sim/Compass_files/mesh_2k/sgr'

temp = 'rd'
#sgr = ["R", "A", "A", "A", "A"]
sgr = ["R", "B", "B", "B", "B"]
hloc = ["112", "112", "132", "122", "142"]

#Find indices of cells along grounding line
dr = xarray.open_dataset(f'{sgr_base}/sgr_data_GLline_x25.nc')
dr.load()
s = np.squeeze(dr.subglacialRunoffFlux.data)
xCell = np.squeeze(dr.xCell.data)
yCell = np.squeeze(dr.yCell.data)

igl = np.where(s > 0)[0]
x = xCell[igl]
y = yCell[igl]
ind = np.argsort(y)
x = x[ind]
y = y[ind]
igl = igl[ind]
print(igl)

xh = x + 20000
iglh = np.zeros(len(x))
x_hor = np.zeros(len(x))
y_hor = np.zeros(len(x))
#--Distance-----------------------
for j in range(len(xh)):
    dist2point = numpy.absolute(numpy.sqrt((xCell-xh[j])**2 + (yCell-y[j])**2))
    ii = numpy.where(dist2point == dist2point.min())[0]
    iglh[j] = np.array(ii)
iglh = iglh.astype(int)
x_hor = xCell[iglh]
y_hor = yCell[iglh]

rho_z = np.zeros((len(sgr), len(igl)))
rho_x = np.zeros((len(sgr), len(igl)))
rho_xz = np.zeros((len(sgr), len(igl)))
rho_zz = np.zeros((len(sgr), len(igl)))
for c in range(len(sgr)):
    #Load files
    fdir = f'{p_base}/{temp}/{temp}_{hloc[c]}{sgr[c]}'

    dsMesh = xarray.open_dataset(f'{fdir}/restart.0003-01-01_00.00.00.nc')
    dsMesh.load()
    minLevelCell = np.squeeze(dsMesh.minLevelCell.data)
    maxLevelCell = np.squeeze(dsMesh.maxLevelCell.data)

    ds = xarray.open_dataset(f'{fdir}/timeSeriesStatsMonthly.0002-12-01.nc')
    ds.load()
    rho = np.squeeze(ds.timeMonthly_avg_potentialDensity.data)

    #Find vertical density gradient along grounding line
    for j in range(len(igl)):
        #print(igl[j])
        #print([j])
        #print(maxLevelCell[igl[j]])
        #print(maxLevelCell[j])
        rho_z[c][j] = (rho[igl[j], minLevelCell[igl[j]]-1] - rho[igl[j], maxLevelCell[igl[j]]-1])
        rho_x[c][j] = (rho[iglh[j], maxLevelCell[iglh[j]] - 1] - rho[igl[j], maxLevelCell[igl[j]] - 1])
        rho_xz[c][j] = (rho[iglh[j], minLevelCell[iglh[j]] - 1] - rho[igl[j], maxLevelCell[igl[j]] - 1])
        rho_zz[c][j] = (rho[iglh[j], minLevelCell[iglh[j]] - 1] - rho[igl[j], minLevelCell[igl[j]] - 1])
        #print(rho_z[c][j])
        #wait = input("Press Enter to continue.")
        #print(rho[:][maxLevelCell-1])

        #rho_x[c][j] = (rho[iglh[j], minLevelCell[j]] - rho[iglh[j], minLevelCell[j]])
        #rho_z[c][j] = (rho[igl[j], minLevelCell[j]] - rho[iglh[j], maxLevelCell[j]])

mkm = 1000
plt.figure(figsize=(4, 4))
for c in range(len(sgr)):
    plt.plot(y/mkm, rho_z[c][:], '-', fillstyle = 'none' , label = f'{hloc[c]}{sgr[c]}')
plt.xlabel('y (km)')
plt.ylabel('rho_z')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()

plt.figure(figsize=(4, 4))
for c in range(len(sgr)):
    plt.plot(y/mkm, rho_x[c][:], '-', fillstyle = 'none' , label = f'{hloc[c]}{sgr[c]}')
plt.xlabel('y (km)')
plt.ylabel('rho_x')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()

plt.figure(figsize=(4, 4))
for c in range(len(sgr)):
    plt.plot(y/mkm, rho_xz[c][:], '-', fillstyle = 'none' , label = f'{hloc[c]}{sgr[c]}')
plt.xlabel('y (km)')
plt.ylabel('rho_xz')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()

plt.figure(figsize=(4, 4))
for c in range(len(sgr)):
    plt.plot(y/mkm, rho_zz[c][:], '-', fillstyle = 'none' , label = f'{hloc[c]}{sgr[c]}')
plt.xlabel('y (km)')
plt.ylabel('rho_zz')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()


plt.figure(figsize=(4, 4))
for c in range(1, len(sgr), 1):
    plt.plot(y/mkm, rho_z[c][:]-rho_z[0][:], '-', fillstyle = 'none' , label = f'{hloc[c]}{sgr[c]}')
plt.xlabel('y (km)')
plt.ylabel('rho_z')
plt.legend(loc = 2)
plt.grid()
plt.rcParams.update({'font.size': 8})
#plt.savefig(f'{dir_fig_save}{site_names[k]}_compare.png', bbox_inches='tight', dpi=300)
plt.show()

for c in range(1, len(sgr), 1):
    print(np.sum(np.squeeze(rho_z[c][:]-rho_z[0][:])))







