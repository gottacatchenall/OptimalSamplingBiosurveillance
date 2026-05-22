import matplotlib.pyplot as plt
from matplotlib.path import Path

import numpy as np
import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM
import cartopy.io.img_tiles as cimgt


def create_map(bbox, zoom):
    quad = cimgt.OSM()
    fig, ax = plt.subplots(
        figsize=(10, 10), 
        subplot_kw={'projection': quad.crs}
    )
    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    ax.add_image(quad, zoom, interpolation='bilinear')
    plt.tight_layout()
    return fig, ax


# circle (bon idx = 10)
bbox = (126.6623598845, 126.67069321780002, 36.77069393740001, 36.77902727070001)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_circle_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()

# triangle (bon idx = 32)
bbox =  (128.52069321040003, 128.5290265437, 36.337360605800015, 36.34569393910001)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_triangle_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()

# x_zoom (bon idx = 29)
bbox =  (127.37069321500002, 127.3790265483, 36.81236060390001, 36.82069393720001)
fig, ax = create_map(bbox, 17)
plt.savefig('plots/korea_x_zoom.svg', dpi=300, bbox_inches='tight')
plt.show()


    