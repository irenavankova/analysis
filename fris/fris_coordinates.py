#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sites_config = [
    {"name": "R01", "lat": -75 - 53.4683 / 60, "lon": -59 - 8.8143 / 60, "type": "apres"},
    {"name": "R02", "lat": -76 - 22.2710 / 60, "lon": -61 - 51.9676 / 60, "type": "apres"},
    {"name": "R03", "lat": -77 - 2.2260 / 60, "lon": -66 - 56.4240 / 60, "type": "apres"},
    {"name": "R04", "lat": -77 - 4.9440 / 60, "lon": -73 - 32.3100 / 60, "type": "apres"},
    {"name": "R05", "lat": -78 - 51.8940 / 60, "lon": -71 - 20.0160 / 60, "type": "apres"},
    {"name": "R06", "lat": -79 - 56.1240 / 60, "lon": -72 - 8.1300 / 60, "type": "apres"},
    {"name": "R07", "lat": -80 - 41.4540 / 60, "lon": -73 - 31.7280 / 60, "type": "apres"},
    {"name": "R08", "lat": -80 - 49.7220 / 60, "lon": -65 - 17.9100 / 60, "type": "apres"},
    {"name": "R09", "lat": -82 - 49.3380 / 60, "lon": -60 - 56.2500 / 60, "type": "apres"},
    {"name": "R10", "lat": -82 - 51.4440 / 60, "lon": -58 - 22.2300 / 60, "type": "apres"},
    {"name": "R12", "lat": -79 - 0.1500 / 60, "lon": -38 - 25.3860 / 60, "type": "apres"},
    {"name": "R13", "lat": -78 - 29.9520 / 60, "lon": -37 - 27.6720 / 60, "type": "apres"},
    {"name": "R14", "lat": -78 - 44.6820 / 60, "lon": -51 - 6.7140 / 60, "type": "apres"},
    {"name": "R15", "lat": -76 - 59.4600 / 60, "lon": -58 - 56.5200 / 60, "type": "apres"},
    {"name": "FSW2", "lat": -80 - 28.8732 / 60, "lon": -44 - 11.3166 / 60, "type": "apres"},
    {"name": "FSE1", "lat": -80 - 58.3560 / 60, "lon": -41 - 26.6640 / 60, "type": "apres"},
    {"name": "FNE1", "lat": -78.54463, "lon": -37.26111, "type": "apres"},
    {"name": "FNE3", "lat": -78.55305, "lon": -38.81890, "type": "apres"},
    {"name": "Site5a", "lat": -80.33333334, "lon": -54.77333305, "type": "apres"},
    {"name": "Site5c", "lat": -80.28672, "lon": -54.70666626, "type": "apres"},
    {"name": "Site2", "lat": -(76 + 42 / 60), "lon": -(64 + 53 / 60), "type": "mooring"},
    {"name": "Site3", "lat": -78.8663, "lon": -71.33753, "type": "mooring"},
    {"name": "Site5", "lat": -(80 + 20 / 60), "lon": -(54 + 46.4 / 60), "type": "mooring"},
    {"name": "Fox1", "lat": -75.5, "lon": -62.908, "type": "fox"},
    {"name": "Fox2", "lat": -75.725, "lon": -61.0557, "type": "fox"},
    {"name": "Fox3", "lat": -75.925, "lon": -59.3, "type": "fox"},
    {"name": "Fox4", "lat": -75.975, "lon": -58.773, "type": "fox"},
    {"name": "FR5", "lat": -75 - 9.8 / 60, "lon": -58 - 43.6 / 60, "type": "fr"},
    {"name": "FR6", "lat": -74 - 42.3 / 60, "lon": -60 - 48.6 / 60, "type": "fr"},
    {"name": "FSW1", "lat": -80.43531, "lon": -44.43134, "type": "mooring"}
]


def find_nearest_mpas_cells(sites, ds_mesh):
    """
    Computes the nearest MPAS grid cell index for a list of coordinate dicts.
    """
    xCell = ds_mesh['xCell'].values
    yCell = ds_mesh['yCell'].values
    zCell = ds_mesh['zCell'].values

    erad_e3 = 6.37122e6
    deg2rad = np.pi / 180.0

    site_names = []
    site_indices = []

    for site in sites:
        # Convert to 0-360 longitude matching MPAS convention
        lon_deg = site["lon"] + 360.0 if site["lon"] < 0 else site["lon"]
        lat_deg = site["lat"]

        site_z = erad_e3 * np.sin(lat_deg * deg2rad)
        site_x = erad_e3 * np.cos(lat_deg * deg2rad) * np.cos(lon_deg * deg2rad)
        site_y = erad_e3 * np.cos(lat_deg * deg2rad) * np.sin(lon_deg * deg2rad)

        dist = np.sqrt((xCell - site_x) ** 2 + (yCell - site_y) ** 2 + (zCell - site_z) ** 2)
        nearest_idx = np.argmin(dist)

        site_names.append(site["name"])
        site_indices.append(nearest_idx)

    return xr.DataArray(site_indices, dims=['site'], coords={'site': site_names})


def plot_mesh_and_sites(sites, ds_mesh):
    """
    Plots the underlying MPAS mesh regional mask (gray for floating FRIS, yellow for Ronne Shelf),
    finds the nearest MPAS cells for the specified sites, and plots them with custom markers.
    """
    # 1. Extract coordinates and transform to degrees
    lat = np.squeeze(ds_mesh.latCell.data) * 180 / np.pi
    lon = np.squeeze(ds_mesh.lonCell.data) * 180 / np.pi
    FloatingMask = np.squeeze(ds_mesh.landIceFloatingMask.data)

    # Calculate Grounding Line mask (iGL) to filter boundaries as done in gmask4TS
    cellsOnCell = ds_mesh.cellsOnCell.values
    nEdgesOnCell = ds_mesh.nEdgesOnCell.values
    nonzero_counts = np.count_nonzero(cellsOnCell, axis=1)
    iGL = (nonzero_counts < nEdgesOnCell)

    # 2. Define masks based on gmask4TS rules
    # FRIS Floating Mask (to be plotted in lightgray)
    iam1 = (FloatingMask == 1) & (lat > -84) & (lat < -74.4) & (lon > 275) & (lon < 331)
    iam2 = (lon < 324.1) | (lat < -78.253)
    fris_floating = np.logical_and(iam1, iam2) & ~iGL

    # 3. Plot background mesh
    plt.figure(figsize=(8, 6))

    # 4. Map site configurations to their closest MPAS elements
    nearest_da = find_nearest_mpas_cells(sites, ds_mesh)

    styles = {
        "apres": {"marker": "s", "color": "blue"},
        "mooring": {"marker": "o", "color": "red"},
        "fox": {"marker": "^", "color": "orange"},
        "fr": {"marker": "D", "color": "green"},
        "none": {"marker": "x", "color": "grey"}
    }

    # Track labels for a clean legend output
    plotted_types = set()

    # 5. Plot nearest cells
    for site in sites:
        style = styles[site["type"]]
        idx = nearest_da.sel(site=site["name"]).values

        # Pull the actual coordinate of the mapped MPAS cell mesh point
        mpas_lon = lon[idx]
        mpas_lat = lat[idx]

        lbl = site["type"] if site["type"] not in plotted_types else ""
        plotted_types.add(site["type"])

        plt.plot(mpas_lon, mpas_lat, marker=style["marker"], color=style["color"],
                 markersize=8, markeredgecolor='black', label=lbl)

    plt.plot(lon[fris_floating], lat[fris_floating], '.', color='black', markersize=2, label='FRIS Floating')

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('MPAS Mesh with Nearest Site Positions')
    plt.legend(loc='upper right', markerscale=1.5)
    plt.grid(True, linestyle='--', alpha=0.5)


if __name__ == "__main__":
    # Path configuration matching your system setup
    fnum = 8
    fris_loc = '/Users/ivankova/Desktop/Fris_hr'
    mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'

    # Open dataset and load required metrics
    dsMesh = xr.open_dataset(mesh_file)
    dsMesh = dsMesh[['latCell', 'lonCell', 'xCell', 'yCell', 'zCell',
                     'landIceFloatingMask', 'cellsOnCell', 'nEdgesOnCell']]
    dsMesh.load()

    # Execute custom mesh mapping and point visualization
    plot_mesh_and_sites(sites_config, dsMesh)
    plt.show()

'''

# Test usage
if __name__ == "__main__":

    # Extract all mooring sites
    moorings = [s for s in sites_config if s["type"] == "mooring"]

    # Extract all Apres sites
    apres = [s for s in sites_config if s["type"] == "apres"]

    # Map types to specific plotting styles
    styles = {
        "apres":   {"marker": "s", "color": "blue"},
        "mooring": {"marker": "o", "color": "red"},
        "fox":     {"marker": "^", "color": "orange"},
        "fr":      {"marker": "D", "color": "green"},
        "none":    {"marker": "x", "color": "grey"}
    }

    # Plot everything cleanly in one go
    for site in sites_config:
        style = styles[site["type"]]
        plt.plot(site["lon"], site["lat"], marker=style["marker"], color=style["color"])

    plt.show()

    all_possible_sites = [
        "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10",
        "R12", "R13", "R14", "R15", "FSW2", "FSE1", "FNE1", "FNE3", "Site5a",
        "Site5c", "Site2", "Site3", "Site5", "Fox1", "Fox2", "Fox3", "Fox4",
        "FR5", "FR6", "FSW1"
    ]
'''


