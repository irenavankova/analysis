#!/usr/bin/env python3

import os
import sys
import numpy as np
import xarray as xr

# Import your custom masking function directly
import gmask_reg


def print_min_max_and_extremes(file_path, mesh_file, variable_name, opt_noGL=1):

    rho_fw = 1000.0
    sec_per_day = 86400.0
    sec_per_year = sec_per_day * 365.0

    try:
        # 1. Use your imported function to generate the mask matrix
        print("Generating FRIS region mask via gmask_reg...")
        iceshelves = ["FRIS"]
        iam = gmask_reg.get_mask(
            iceshelves, mesh_file, opt_noGL=opt_noGL, opt_wct=0
        )
        fris_mask = iam[0, :]  # Extract the FRIS boolean array

        # 2. Open the primary data NetCDF file
        with xr.open_dataset(file_path) as ds:

            if variable_name not in ds:
                print(f"Error: Variable '{variable_name}' not found in the file.")
                print(f"Available variables: {list(ds.data_vars)}")
                return

            # Extract data array
            data_array = ds[variable_name]

            if variable_name == "timeMonthly_avg_landIceFreshwaterFluxTotal":
                data_array = data_array * sec_per_year / rho_fw

            # Convert to numpy array
            data_np = data_array.values

            # Apply the spatial FRIS mask
            # If data_np has a Time dimension (e.g., shape: [Time, nCells]), mask along axis 1
            if data_np.ndim == 2:
                masked_data = data_np[:, fris_mask]
            else:
                masked_data = data_np[fris_mask]

            # Flatten and drop NaNs/missing values
            data_flat = masked_data.flatten()
            data_flat = data_flat[~np.isnan(data_flat)]

            if len(data_flat) == 0:
                print("Error: No valid numeric data found within the FRIS mask.")
                return

            # Determine how many elements to grab (handles cases with < 20 points)
            num_elements = min(50, len(data_flat))

            # Efficiently pull and sort the lowest values
            lowest_20 = np.partition(data_flat, num_elements)[:num_elements]
            lowest_20.sort()

            # Efficiently pull and sort the highest values
            highest_20 = np.partition(data_flat, -num_elements)[-num_elements:]
            highest_20.sort()

            # Print results to the screen
            gl_status = "WITHOUT Grounding Line" if opt_noGL else "WITH Grounding Line"
            print("=" * 60)
            print(f"Variable: {variable_name}")
            print(f"Region:   FRIS ({gl_status})")
            print(f"Total Grid Cells Analyzed: {len(data_flat)}")
            print("=" * 60)
            print(f"Absolute Minimum Value: {lowest_20[0]}")
            print(f"Absolute Maximum Value: {highest_20[-1]}")
            print("-" * 60)

            print(f"--- LOWEST {num_elements} VALUES (Ascending) ---")
            for i, val in enumerate(lowest_20, 1):
                print(f"{i:2d}: {val}")

            print(f"\n--- HIGHEST {num_elements} VALUES (Descending) ---")
            for i, val in enumerate(reversed(highest_20), 1):
                print(f"{i:2d}: {val}")
            print("=" * 60)

    except FileNotFoundError as e:
        print(f"Error: A required file could not be found.\n{e}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    Fnum = "1"

    # Paths configured from your examples
    FILE_NAME = f"/Users/ivankova/Desktop/Fris_hr/Fris_ncfiles/F{Fnum}/ncfiles/F{Fnum}melt_annual_2D_Y4.nc"
    MESH_FILE = f"/Users/ivankova/Desktop/Fris_hr/Fris_ncfiles/F{Fnum}/ncfiles/F{Fnum}mesh.nc"

    VARIABLE = "timeMonthly_avg_landIceFreshwaterFluxTotal"

    # Set opt_noGL=1 to exclude Grounding Line cells, or opt_noGL=0 to include them
    print_min_max_and_extremes(
        FILE_NAME, MESH_FILE, VARIABLE, opt_noGL=1
    )