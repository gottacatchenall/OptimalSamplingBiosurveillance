import matplotlib.pyplot as plt
from matplotlib.path import Path

import numpy as np
import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM
import cartopy.io.img_tiles as cimgt


def create_map(bbox, zoom):
    quad = cimgt.QuadtreeTiles()
    fig, ax = plt.subplots(
        figsize=(10, 10), 
        subplot_kw={'projection': quad.crs}
    )
    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    ax.add_image(quad, zoom, interpolation='bilinear')
    plt.tight_layout()
    return fig, ax

# circle (bon idx = 3)
bbox = (126.8206932172, 126.8290265505, 36.4373606054, 36.4456939387)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_circle_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()

# triangle (bon idx = 5)
bbox =  (128.05402654559998, 128.0623598789, 35.8623606077, 35.870693941)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_triangle_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()

# x_zoom (bon idx = 15)
bbox =  (126.5706932182, 126.5790265515, 36.904027270200004, 36.912360603500005)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_x_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()


