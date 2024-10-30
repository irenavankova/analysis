#!/usr/bin/env python3
import numpy as np
import xarray
import numpy # for arrays!
#import scipy.io  # load matfiles
import matplotlib.pyplot as plt
import gsw
import os

def get_mask(is_list):

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

    #H = np.squeeze(dsMesh.restingThickness.data)
    #H = np.nansum(H, axis=1); H = np.squeeze(H)

    iam = np.zeros((len(is_list), len(FloatingMask)))

    for n in range(len(is_list)):
        if is_list[n] == "Amery":
            #AMERY
            iam[n,:] = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
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
        elif is_list[n] == "Ekstrom":
            # Ekstrom
            iam[n,:] = (FloatingMask == 1) & (lat > -72) & (lat < -70) & (lon > 350) & (lon < 352.4)
        elif is_list[n] == "Thwaites":
            # Thwaites
            iam[n,:] = (FloatingMask == 1) & (lat > -75.6894) & (lat < -74.6915) & (lon > -107.9729 + 360) & (lon < -104.0573 + 360)
        elif is_list[n] == "Pine Island":
            # Pine Island
            iam[n,:] = (FloatingMask == 1) & (lat > -75.7071) & (lat < -74.2684) & (lon > -103.4880 + 360) & (lon < -98.4832 + 360)
        elif is_list[n] == "Getz":
            # Getz
            iam[n, :] = (FloatingMask == 1) & (lat > -75) & (lat < -73.5) & (lon > 225) & (lon < 245.5)


    #iam = iam.astype(int)
    iam = iam.astype(bool)
    return iam, areaCell