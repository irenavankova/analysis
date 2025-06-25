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
        if is_list[n] == "Antarctica":
            #Antarctica
            iam[n,:] = (FloatingMask == 1)
        elif is_list[n] == "Belli":
            # Peninnsula West -Bellingshausen
            iam1 = (FloatingMask == 1) & (lat > -77.4) & (lon > 261) & (lon < 293.8)
            iam2 = (FloatingMask == 1) & (lat > -77.3) & (lon > 255) & (lon < 262.8)
            iam[n,:] = np.logical_or(iam1, iam2)
        elif is_list[n] == "Amundsen":
            # Amundsen (up to Getz)
            iam[n,:] = (FloatingMask == 1) & (lat > -76) & (lat < -73.2) & (lon > 225) & (lon < 261)
        elif is_list[n] == "Ross":
            # Ross
            iam1 = (FloatingMask == 1) & (lat > -85.6) & (lat < -77.4) & (lon > 158.64) & (lon < 212.5)
            iam2 = (lat < -77.8) | (lon < 200)
            iam[n,:] = np.logical_and(iam1, iam2)
        elif is_list[n] == "Eastant":
            # East Antarctica up to Amery
            iam[n,:] = (FloatingMask == 1) & (lat > -72) & (lon > 80) & (lon < 166)
        elif is_list[n] == "Amery":
            # Amery
            iam[n,:] = (FloatingMask == 1) & (lat > -74.3608) & (lat < -67.7122) & (lon > 62.0419) & (lon < 78)
        elif is_list[n] == "Dml":
            # DML
            iam[n,:] = (FloatingMask == 1) & ((lon < 62.) | (lon > 332))
        elif is_list[n] == "Fris":
            # FRIS
            iam[n,:] = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)
        elif is_list[n] == "Larsens":
            # Peninnsula East -Larsens
            iam[n,:] = (FloatingMask == 1) & (lat > -74.3) & (lon > 294) & (lon < 301)
        elif is_list[n] == "si_so":
            # Southern Ocean
            iam[n,:] = (lat < -50) & (FloatingMask == 0)
        elif is_list[n] == "si_ea":
            # East Ant
            iam[n,:] = (lat < -60) & (lon > 90) & (lon < 150) & (FloatingMask == 0)
        elif is_list[n] == "si_amery":
            # Amery
            iam[n,:] = (lat < -60) & (lon > 60) & (lon < 90) & (FloatingMask == 0)
        elif is_list[n] == "si_dml":
            # DML
            iam1 = (lat < -60) & (lon > 350) & (FloatingMask == 0)
            iam2 = (lat < -60) & (lon > 0) & (lon < 60) & (FloatingMask == 0)
            iam[n,:] = np.logical_or(iam1, iam2)
        elif is_list[n] == "si_weddell":
            # Weddell
            iam[n,:] = (lat < -60) & (lon > 297) & (lon < 350) & (FloatingMask == 0)
        elif is_list[n] == "si_ab":
            # AB
            iam[n,:] = (lat < -68) & (lon > 220) & (lon < 294) & (FloatingMask == 0)
        elif is_list[n] == "si_rosse":
            # Ross East
            iam[n,:] = (lat < -60) & (lon > 180) & (lon < 220) & (FloatingMask == 0)
        elif is_list[n] == "si_rossw":
            # Ross West
            iam[n,:] = (lat < -60) & (lon > 150) & (lon < 180) & (FloatingMask == 0)

    iam = iam.astype(bool)
    return iam