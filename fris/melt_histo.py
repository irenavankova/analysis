#!/usr/bin/env python3

import xarray as xr
import numpy as np
import gmask_reg
import os

import matplotlib.pyplot as plt


fy = '2-4'
opt_plot = 2
opt_save = 1
fris_loc = '/Users/irenavankova/Desktop/Fris_hr'

rho_fw = 1000.
sec_per_day = 86400.
sec_per_year = sec_per_day * 365.

# Define the list of values to loop over
fnums = [1, 2, 4, 8]
savename = 'F1248'
iceshelves = ["Fris","Evans","Institute","Rutford","Foundation","Support Force","Recovery","Slessor"]
#iceshelves = ["Rutford"]

fig, axes = plt.subplots(3, 3, figsize=(16, 12))
axes = axes.flatten() # Flatten to 1D array for easy indexing

for i, shelf_name in enumerate(iceshelves):
    ax = axes[i]
    for fnum in fnums:
        mesh_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}mesh.nc'
        out_file = f'{fris_loc}/Fris_ncfiles/F{fnum}/ncfiles/F{fnum}melt_annual_2D_Y{fy}.nc'

        # Get mask and indices (assuming gmask_reg is already defined in your environment)
        iam = gmask_reg.get_mask([shelf_name], mesh_file)
        iis = iam[0, :]

        # Load Mesh and select indices
        dsMesh = xr.open_dataset(mesh_file).isel(nCells=iis)
        FloatingMask = np.squeeze(dsMesh['landIceFloatingMask'].values)
        areaCell = np.squeeze(dsMesh['areaCell'].values)

        # Load Data and select indices
        ds = xr.open_dataset(out_file).isel(nCells=iis)
        lifw = ds['timeMonthly_avg_landIceFreshwaterFluxTotal'].values.flatten()

        # Calculate Melt rate (m/yr)
        shelf_indices = np.where(FloatingMask == 1)[0]
        area_shelf = areaCell[shelf_indices]
        melt_shelf = lifw[shelf_indices] * sec_per_year / rho_fw

        total_shelf_area = np.sum(area_shelf)
        melt_weighted_mean = np.sum(melt_shelf * area_shelf) / total_shelf_area

        plot_range = (-7.5, 7.5)
        n_bins = 60
        bin_width = (plot_range[1] - plot_range[0]) / n_bins

        # 4. Calculate weights for the shelf cells
        weights = area_shelf

        ax.hist(melt_shelf, bins=n_bins, range=plot_range, weights=weights,
                 alpha=0.8, label=f'F{fnum} (Mean melt: {melt_weighted_mean:.2f} m/yr)',
                 histtype='step', linewidth=2)

    ax.set_title(f'{shelf_name}', fontsize=12)
    if i >= 5: ax.set_xlabel('Melt Rate (m/yr)')
    if i % 3 == 0: ax.set_ylabel('Melt * cell Area')
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)

for j in range(len(iceshelves), len(axes)):
    axes[j].axis('off')

plt.tight_layout()

if opt_save == 1:
    plt.savefig(f'{fris_loc}/Fris_plots/melt/Melt_hist_{savename}_Y{fy}.png',
                bbox_inches='tight', dpi=300)
else:
    plt.show()