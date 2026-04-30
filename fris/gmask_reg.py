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
        if is_list[n] == "FRIS":
            # FRIS
            #iam[n,:] = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)

            iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
            iam2 = (lon < 324.1) | (lat < -78.253)
            iam[n,:] = np.logical_and(iam1, iam2)
        elif is_list[n] == "Evans":
            # FRIS Evans
            iam[n,:] = (FloatingMask == 1) & (lat > -77.2) & (lat < -75.75) & (lon > 281) & (lon < 285)
        elif is_list[n] == "Rutford":
            # FRIS Rutford
            iam[n,:] = (FloatingMask == 1) & (lat > -79.6) & (lat < -78.2) & (lon > 275.5) & (lon < 281)
        elif is_list[n] == "Institute":
            # FRIS Institute
            iam[n,:] = (FloatingMask == 1) & (lat > -81) & (lat < -80.5) & (lon > 284) & (lon < 289)
        elif is_list[n] == "Foundation":
            # FRIS Foundation
            iam[n,:] = (FloatingMask == 1) & (lat > -84) & (lat < -83.04) & (lon > 298) & (lon < 305)
        elif is_list[n] == "Support Force":
            # FRIS Support Force
            iam[n,:] = (FloatingMask == 1) & (lat > -83) & (lat < -81.8) & (lon > 313.5) & (lon < 317)
        elif is_list[n] == "Recovery":
            # FRIS Recovery
            iam[n,:] = (FloatingMask == 1) & (lat > -81) & (lat < -80.62) & (lon > 321) & (lon < 325)
        elif is_list[n] == "Slessor":
            # FRIS Slessor
            iam[n,:] = (FloatingMask == 1) & (lat > -80.5) & (lat < -79.9) & (lon > 327) & (lon < 332)
        elif is_list[n] == "Bailey":
            # FRIS Bailey
            iam[n, :] = (FloatingMask == 1) & (lat > -79.7) & (lat < -79) & (lon > 326) & (lon < 332)

    iam = iam.astype(bool)
    return iam