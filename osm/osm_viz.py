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




# Okay we're getting sentinel 2 with the planetary computer and overlaying OSM. 

import pystac_client
import planetary_computer
import rasterio
import rioxarray
import xarray
import stackstac

# Open the STAC API
catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace
)

bbox = (126.8206932172, 36.4373606054, 126.8290265505, 36.4456939387)

# Search for Sentinel-2 Level-2A data
search = catalog.search(
    collections=["sentinel-2-l2a"],
    bbox=bbox,
    datetime="2023-06-01/2023-09-30",
    query={"eo:cloud_cover": {"lt": 10}}
)

first_item = next(search.items())


def read_band(asset, band_id, crs = 4326):
    with rasterio.open(asset[band_id].href) as src:
        with WarpedVRT(src, crs=crs) as vrt:
            dst_window = vrt.window(*bbox)
            data = vrt.read(window=dst_window)
            return data
        

stack = stackstac.stack(first_item, epsg=4326).sel(band=["B04", "B03", "B02"])

xr_stack = stackstac.stack(first_item, epsg=4326, bounds_latlon=bbox)


a = xr_stack.sel(band="B04").compute()

    