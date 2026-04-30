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
kgInGt = 1e12

xmin = -3.5
xmax = 7.5
n_bins = 20

# Define the list of values to loop over
fnums = [8,4,2,1]
savename = 'F1248'
iceshelves = ["FRIS","Evans","Institute","Rutford","Foundation","Support Force","Recovery","Slessor","Bailey"]
#iceshelves = ["Rutford"]
color_map = {
    1: 'black',
    2: 'dodgerblue',
    4: 'orange',
    8: 'brown'
}

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
        total_melt_flux = np.sum(melt_shelf * area_shelf)*rho_fw/kgInGt

        plot_range = (xmin, xmax)
        #n_bins = 60
        #bin_width = (plot_range[1] - plot_range[0]) / n_bins
        bins = np.linspace(plot_range[0], plot_range[1], n_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_width = bins[1] - bins[0]

        # 4. Calculate weights for the shelf cells
        weights = area_shelf

        # 3. Compute histogram values (counts) manually
        # This returns the heights for each bin
        hist_vals, _ = np.histogram(melt_shelf, bins=bins, weights=weights)

        # 4. Plot as a line connecting bin centers
        current_color = color_map.get(fnum, 'gray')
        ax.plot(bin_centers, hist_vals,
                label=f'F{fnum} (Melt flux: {total_melt_flux:.2f} Gt/yr)',
                color=current_color,
                linewidth=1.5,
                alpha=1.0)

        # Optional: Fill the area under the line slightly for better visibility
        #ax.fill_between(bin_centers, 0, hist_vals, alpha=0.1)
        ax.axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.set_xlim([xmin, xmax])

    ax.set_title(f'{shelf_name}', fontsize=12)
    if i >= 5: ax.set_xlabel('Melt Rate (m/yr)')
    if i % 3 == 0: ax.set_ylabel('Volume flux through cell')
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)

for j in range(len(iceshelves), len(axes)):
    axes[j].axis('off')

plt.tight_layout()

if opt_save == 1:
    plt.savefig(f'{fris_loc}/Fris_plots/melt/Melt_hist_V2_{savename}_Y{fy}.png',
                bbox_inches='tight', dpi=300)
else:
    plt.show()