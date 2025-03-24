#!/usr/bin/env python3
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Create a polar stereographic projection for Antarctica
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Stereographic(central_longitude=0, central_latitude=-90))

# Set the extent of the map to focus on the southern polar region
ax.set_extent([-180, 180, -90, 0], crs=ccrs.PlateCarree())

# Add coastlines and countries (this will be drawn only for the land in the Antarctic region)
#ax.coastlines(resolution='110m', color='black', linewidth=1)
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightblue')
ax.add_feature(cfeature.OCEAN, edgecolor='none')

# Optionally, add gridlines for better visualization
ax.gridlines(draw_labels=True)

# Title and labels
ax.set_title('Masked Antarctic Continent on a Polar Stereographic Projection')

# Display the plot
plt.show()