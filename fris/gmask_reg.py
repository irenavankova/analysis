#!/usr/bin/env python3
import numpy as np
import xarray

def get_mask(is_list,ocean_rst_file):
    #is_list is list of region names
    #ocean_rst_file is MPAS-ocean file

    dsMesh = xarray.open_dataset(ocean_rst_file)
    dsMesh = dsMesh[['latCell', 'lonCell','landIceFloatingMask']]
    dsMesh.load()

    lat = np.squeeze(dsMesh.latCell.data)
    lon = np.squeeze(dsMesh.lonCell.data)
    lat = lat*180/np.pi
    lon = lon*180/np.pi
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

    iam = np.zeros((len(is_list), len(FloatingMask)))

    for n in range(len(is_list)):
        if is_list[n] == "Fris":
            # FRIS
            #iam[n,:] = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)

            iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
            iam2 = (lon < 324.1) | (lat < -78.253)
            iam[n,:] = np.logical_and(iam1, iam2)

    iam = iam.astype(bool)
    return iam