# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 21:30:22 2021

@author: sbe002
This is to test gridgeo as in 
https://nbviewer.org/github/pyoceans/gridgeo/blob/master/notebooks/gridgeo_tour.ipynb
source: https://github.com/pyoceans/gridgeo
"""
import gridgeo
import netCDF4
import os
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist

dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
file_nc_input = os.path.join(dir_testinput,'input','westerscheldt_net.nc')
file_nc_map = os.path.join(dir_testinput,'output','westerscheldt01_map.nc')
file_nc_his = os.path.join(dir_testinput,'output','westerscheldt01_his.nc')

ds = netCDF4.Dataset(file_nc_map)
#get lists with vars/dims, times, station/crs/structures
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_input)

grid = gridgeo.GridGeo(file_nc_input, mesh='mesh2d')


import fiona

schema = {
    'geometry': 'MultiPolygon',
    'properties': {'name': f'str:{len(grid.mesh)}'}
}

with fiona.open('grid.shp', 'w', 'ESRI Shapefile', schema) as f:
    f.write(
        {
            'geometry': grid.__geo_interface__,
            'properties': {'name': grid.mesh}
        }
    )


[s for s in dir(grid) if not s.startswith('_')]

grid.mesh

print(f'The grid has {len(grid.geometry)} polygons, showing the first 5.')

grid.geometry[:5]

grid.outline

import cartopy.crs as ccrs
import matplotlib.pyplot as plt


fig, ax = plt.subplots(
    figsize=(12, 12),
    subplot_kw={'projection': ccrs.PlateCarree()}
)

kw = dict(linestyle='-', alpha=0.25, color='darkgray')
ax.triplot(grid.triang, **kw)
ax.coastlines(resolution='10m');

kw = {
    'fill': '#fd7d11',
    'fill_opacity': 0.2,
    'stroke_opacity': 1,
    'float_precision': 2,
}

geojson = grid.to_geojson(**kw)
geojson['properties']

grid.__geo_interface__.keys()

## Saving the grid to as geojson file
grid.save('grid.geojson', **kw)

## shapefile
grid.save('grid.shp')