#!/usr/bin/env python3
import numpy as np
import xarray

def get_mask(is_list,ocean_rst_file,opt_noGL=1, opt_wct = 1):
    #is_list is list of region names
    #ocean_rst_file is MPAS-ocean file

    dsMesh = xarray.open_dataset(ocean_rst_file)
    if opt_wct == 1:
        dsMesh = dsMesh[['latCell', 'lonCell', 'landIceFloatingMask', 'cellsOnCell', 'nEdgesOnCell', 'layerThickness']]
    else:
        dsMesh = dsMesh[['latCell', 'lonCell', 'landIceFloatingMask', 'cellsOnCell', 'nEdgesOnCell']]

    dsMesh.load()

    lat = np.squeeze(dsMesh.latCell.data)
    lon = np.squeeze(dsMesh.lonCell.data)
    lat = lat*180/np.pi
    lon = lon*180/np.pi
    FloatingMask = np.squeeze(dsMesh.landIceFloatingMask.data)

    if opt_wct == 1:
        wct = dsMesh['layerThickness'].sum(dim='nVertLevels')
        wct = np.squeeze(wct)


    iam = np.zeros((len(is_list), len(FloatingMask)))

    cellsOnCell = dsMesh.cellsOnCell.values
    cellsOnCell = cellsOnCell

    nEdgesOnCell = dsMesh.nEdgesOnCell.values
    nEdgesOnCell = nEdgesOnCell

    nonzero_counts = np.count_nonzero(cellsOnCell, axis=1)
    # iGL = np.where(nonzero_counts < nEdgesOnCell)[0]
    iGL = (nonzero_counts < nEdgesOnCell)

    for n in range(len(is_list)):
        #current_mask = np.zeros((len(FloatingMask)))
        current_mask = np.zeros((len(FloatingMask)), dtype=bool)
        if is_list[n] == "FRIS":
            # FRIS
            #iam[n,:] = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 330)

            iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
            iam2 = (lon < 324.1) | (lat < -78.253)
            current_mask = (iam1 & iam2)
        elif is_list[n] == "Evans":
            # FRIS Evans
            current_mask = (FloatingMask == 1) & (lat > -77.2) & (lat < -75.75) & (lon > 281) & (lon < 285)
        elif is_list[n] == "Rutford":
            # FRIS Rutford
            current_mask = (FloatingMask == 1) & (lat > -79.6) & (lat < -78.2) & (lon > 275.5) & (lon < 281)
        elif is_list[n] == "Institute":
            # FRIS Institute
            current_mask = (FloatingMask == 1) & (lat > -81) & (lat < -80.5) & (lon > 284) & (lon < 289)
        elif is_list[n] == "Foundation":
            # FRIS Foundation
            current_mask = (FloatingMask == 1) & (lat > -84) & (lat < -83.04) & (lon > 298) & (lon < 305)
        elif is_list[n] == "Support Force":
            # FRIS Support Force
            current_mask = (FloatingMask == 1) & (lat > -83) & (lat < -81.8) & (lon > 313.5) & (lon < 317)
        elif is_list[n] == "Recovery":
            # FRIS Recovery
            current_mask = (FloatingMask == 1) & (lat > -81) & (lat < -80.62) & (lon > 321) & (lon < 325)
        elif is_list[n] == "Slessor":
            # FRIS Slessor
            current_mask = (FloatingMask == 1) & (lat > -80.5) & (lat < -79.9) & (lon > 327) & (lon < 332)
        elif is_list[n] == "Bailey":
            # FRIS Bailey
            current_mask = (FloatingMask == 1) & (lat > -79.7) & (lat < -79) & (lon > 326) & (lon < 332)
        elif is_list[n] == "FRISshelf":
            # FRIS Shelf
            current_mask = (FloatingMask == 0) & (lat > -80) & (lat < -72) & (lon > 298) & (lon < 332) & (wct < 1500)
        elif is_list[n] == "RonneDshelf":
            # Ronne depression shelf
            current_mask = (FloatingMask == 0) & (lat > -80) & (lat < -72.5) & (lon > 298) & (lon < 304.5) & (wct < 1200) & (wct > 500)
        elif is_list[n] == "FilchnerDshelf":
            # Filchner depression shelf
            current_mask = (FloatingMask == 0) & (lat > -80) & (lat < -74.9) & (lon > 316) & (lon < 332) & (wct < 1200) & (wct > 500)
        elif is_list[n] == "BerknerBank":
            # Berkner Bank
            current_mask = (FloatingMask == 0) & (lat > -80) & (lat < -76) & (lon > 307) & (lon < 316) & (wct < 300)
        elif is_list[n] == "RonneDcavity":
            # Ronne depression cavity
            current_mask = (FloatingMask == 1) & (lat > -79) & (lat < -75.5) & (lon > 282) & (lon < 294) & (wct > 500)
        elif is_list[n] == "FilchnerDcavity":
            # Filchner depression cavity
            current_mask = (FloatingMask == 1) & (lat > -81.5) & (lat < -78) & (lon > 317) & (lon < 324) & (wct > 500)
        elif is_list[n] == "BerknerSouth":
            # Berkner South
            current_mask = (FloatingMask == 1) & (lat > -84) & (lat < -80.5) & (lon > 300) & (lon < 316) & (wct > 500)
        elif is_list[n] == "Shelf":
            # WeddellSouth
            #current_mask = (lat > -84) & (lat < -70) & (lon > 275) & (lon < 335) & (wct < 2000)
            current_mask = (wct < 2000)

        if opt_noGL == 1:
            current_mask = current_mask & ~iGL

        iam[n, :] = current_mask

    iam = iam.astype(bool)
    return iam

# Test usage
if __name__ == "__main__":

    iceshelves = ["FRIS", "Rutford"]
    mesh_file = '/Users/ivankova/Desktop/Fris_hr/E3SM_init/ocean.ECwISC30to60E2r1.230220.nc'

    print("-----GL")
    iam = get_mask(iceshelves, mesh_file, opt_noGL=0)
    iis = iam[0, :]
    bool_array = np.array(iis)

    # Print the array and its shape
    print("Shape GL:", bool_array.shape)
    count = np.sum(bool_array)

    print("Number of cells where array == 1:", count)

    print("-----no GL")
    iam = get_mask(iceshelves, mesh_file, opt_noGL=1)
    iis = iam[0, :]
    bool_array = np.array(iis)

    # Print the array and its shape
    print("Shape GL:", bool_array.shape)
    count = np.sum(bool_array)

    print("Number of cells where array == 1:", count)



