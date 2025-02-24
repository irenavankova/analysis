#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

def get_mask(is_list):

    #MPAS-ocean file
    p_file = f'/Users/irenavankova/Work/data_sim/E3SM_files/E3SM_initial_condition/SOwISC12to60E2r4/ocean.SOwISC12to60E2r4.230220.nc'

    dsMesh = xarray.open_dataset(p_file)
    dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask','nVertLevels','areaCell','maxLevelCell','restingThickness']]
    dsMesh.load()

    lat = np.squeeze(dsMesh.latCell.data)
    lon = np.squeeze(dsMesh.lonCell.data)
    lat = lat*180/np.pi
    lon = lon*180/np.pi
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)
    areaCell = np.squeeze(dsMesh.areaCell.data)

    H = np.squeeze(dsMesh.restingThickness.data)
    H = np.nansum(H, axis=1); H = np.squeeze(H)

    iam = np.zeros((len(is_list), len(FloatingMask)))

    # Ship CTD data: Shenjie Zhou
    sz_file = f'/Users/irenavankova/Desktop/Shenjie/Shenjie_CT.nc'
    dsSZ = xarray.open_dataset(sz_file)
    dsSZ = dsSZ[['lat', 'lon', 'bathy']]
    dsSZ.load()
    latsz = np.squeeze(dsSZ.lat.data)
    lonsz = np.squeeze(dsSZ.lon.data)
    lonsz[lonsz < 0] = lonsz[lonsz < 0] + 360
    Hsz = np.squeeze(dsSZ.bathy.data)

    isz = np.zeros((len(is_list), len(latsz)))

    for n in range(len(is_list)):
        if is_list[n] == "Amery":
            #AMERY
            iam[n,:] = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
        elif is_list[n] == "Amery_shelf":
            #AMERY continental shelf
            iam[n,:] = (FloatingMask == 0) & (lat > -70) & (lat < -65) & (lon > 67.5) & (lon < 80) & (H < 1500)
            isz[n,:] = (latsz > -70) & (latsz < -65) & (lonsz > 67.5) & (lonsz < 80) & (-Hsz < 1500)
        elif is_list[n] == "RoiB":
            # Roi Baudouin
            iam[n,:] = (FloatingMask == 1) & (lat > -71.5) & (lat < -69) & (lon > 24) & (lon < 33)
        elif is_list[n] == "Munin":
            # Muninisen
            iam[n,:] = (FloatingMask == 1) & (lat > -71.5) & (lat < -69) & (lon > 19) & (lon < 22)
        elif is_list[n] == "Nivl":
            # Nivlisen
            iam[n,:] = (FloatingMask == 1) & (lat > -71) & (lat < -69.8) & (lon > 9.5) & (lon < 12.9)
        elif is_list[n] == "Fimbul":
            # Fimbul
            iam[n,:] = (FloatingMask == 1) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8))
        elif is_list[n] == "Fimbul_shelf":
            # Fimbul shelf
            iam[n, :] = (FloatingMask == 0) & (lat > -72) & (lat < -69.5) & ((lon > 357.3) | (lon < 7.8)) & (H < 1500)
            isz[n, :] = (latsz > -72) & (latsz < -69.5) & ((lonsz > 357.3) | (lonsz < 7.8)) & (-Hsz < 1500)
        elif is_list[n] == "Ekstrom":
            # Ekstrom
            iam[n,:] = (FloatingMask == 1) & (lat > -72) & (lat < -70) & (lon > 350) & (lon < 352.4)
        elif is_list[n] == "Thwaites":
            # Thwaites
            iam[n,:] = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729 + 360) & (lon < -104.0573 + 360)
        elif is_list[n] == "Pine_Island":
            # Pine Island
            iam[n,:] = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880 + 360) & (lon < -98.4832 + 360)
        elif is_list[n] == "Getz":
            # Getz
            iam[n, :] = (FloatingMask == 1) & (lat > -75) & (lat < -73.5) & (lon > 225) & (lon < 245.5)
        elif is_list[n] == "Totten":
            # Totten
            iam[n, :] = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 117.5)
        elif is_list[n] == "MoscowU":
            # Moscow University
            iam[n, :] = (FloatingMask == 1) & ((lat > -68) & (lat < -66) & (lon > 119.5) & (lon < 122.3)) & ((lat < -67) | (lon > 120.5))
        elif is_list[n] == "TottenMU":
            # Totten MoscowU
            iam[n, :] = (FloatingMask == 1) & (lat > -68) & (lat < -66) & (lon > 113.5) & (lon < 122.3)
        elif is_list[n] == "TottenMU_shelf":
            # Totten MoscowU shelf
            iam[n, :] = (FloatingMask == 0) & (lat > -68) & (lat < -65) & (lon > 114) & (lon < 124) & (H < 1500)
            isz[n, :] = (latsz > -68) & (latsz < -65) & (lonsz > 114) & (lonsz < 124) & (-Hsz < 1500)
        elif is_list[n] == "George_VI":
            # GeorgeVI
            iam[n, :] = (FloatingMask == 1) & ((lat > -74) & (lat < -70) & (lon > 287.3) & (lon < 293.4)) & ((lat < -72.5) | (lon > 290))
        elif is_list[n] == "Stange":
            # Stange
            iam[n, :] = (FloatingMask == 1) & (lat > -73.5) & (lat < -72.5) & (lon > 271) & (lon < 274)
        elif is_list[n] == "Abbot":
            # Abbot
            iam[n, :] = (FloatingMask == 1) & ((lat > -73.34) & (lat < -2) & (lon > 256) & (lon < 271)) & ((lat < -72.28) | (lon < 260)) & ((lat > -73.2) | (lon > 260))
        elif is_list[n] == "Larsen_C":
            # Larsen
            iam[n, :] = (FloatingMask == 1) & (lat > -69.5) & (lat < -66.1) & (lon > 294) & (lon < 300)
        elif is_list[n] == "Larsen_C_shelf":
            # Larsen shelf
            iam[n, :] = (FloatingMask == 0) & (lat > -69.5) & (lat < -66.1) & (lon > 298) & (lon < 310) & (H < 1500)
            isz[n, :] = (latsz > -69.5) & (latsz < -66.1) & (lonsz > 298) & (lonsz < 310) & (-Hsz < 1500)
        elif is_list[n] == "Filchner-Ronne":
            # FRIS
            iam[n, :] = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
        elif is_list[n] == "Filchner-Ronne_shelf":
            # FRIS
            iam[n, :] = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (H < 1500)
            isz[n, :] = (latsz > -80) & (latsz < -72) & (lonsz > 298) & (lonsz < 332) & (-Hsz < 1500)
        elif is_list[n] == "Ross":
            # Ross
            iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
            iam2 = (lat < -77.8) | (lon < 200)
            iam[n, :] = np.logical_and(iam1, iam2)
        elif is_list[n] == "Ross_shelf":
            # Ross shelf
            iam[n, :] = (FloatingMask == 0) & (lat > -85.6) & (lat < -73) & (lon > 158.64) & (lon < 210) & (H < 1500)
            isz[n, :] = (latsz > -85.6) & (latsz < -73) & (lonsz > 158.64) & (lonsz < 210) & (-Hsz < 1500)
        elif is_list[n] == "Amundsen":
            # Amundsen sea ice shelves all
            iam[n, :] = (FloatingMask == 1) & (lat > -76) & (lat < -73.2) & (lon > 225) & (lon < 261)
        elif is_list[n] == "Amundsen_shelf":
            iam[n, :] = (FloatingMask == 0) & (lat > -76) & (lat < -71) & (lon > 225) & (lon < 260) & (H < 1500)
            isz[n, :] = (latsz > -76) & (latsz < -71) & (lonsz > 225) & (lonsz < 260) & (-Hsz < 1500)

    #iam = iam.astype(int)
    iam = iam.astype(bool)
    isz = isz.astype(bool)
    return iam, areaCell, isz