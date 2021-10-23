# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 22:29:12 2021

@author: sbe002
from: https://gis.stackexchange.com/questions/369026/extract-shapely-polygon-from-matplotlib-hexbin-polycollection
"""
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from matplotlib.patches import RegularPolygon


geoms = [Point(919.000, 533.500), Point(940.500, 528.500), 
 Point(937.000, 555.000), Point(889.000, 526.500), Point(898.000, 520.500)]

l = gpd.GeoDataFrame({'geometry':geoms})

y= np.array([ 533.500,  528.500,  555.000,  526.500,  520.500])
x= np.array([ 919.000,  940.500,  937.000,  889.000,  898.000])

f,ax=plt.subplots(1)
im =ax.hexbin(x,y, gridsize=(3,3), edgecolor='black')
l.plot(ax=ax, color='black', markersize=55)

array_of_hexes=[]
for x,y in im.get_offsets():
    hexes = RegularPolygon((x, y), numVertices=6, radius= 5 )
    verts = hexes.get_path().vertices
    trans = im.get_transform()
    points = trans.transform(verts)
    array_of_hexes.append(Polygon(points))

gpd.GeoDataFrame({'geometry':array_of_hexes}).plot(edgecolor='black')


## Incase anyone faces a similar issue, the resolution is:
# collection = plt.hexbin(x,y, gridsize=(5,5))
collection = plt.hexbin(x, y, gridsize=(5,5) )
hex_polys = collection.get_paths()[0].vertices
hex_array = []
for xs,ys in collection.get_offsets():
    hex_x = np.add(hex_polys[:,0],  xs)
    hex_y = np.add(hex_polys[:,1],  ys)
    hex_array.append(Polygon(np.vstack([hex_x, hex_y]).T))

hex_grid = gpd.GeoDataFrame({'geometry':hex_array})

 
# Make the plot
plt.hexbin(x, y, gridsize=(5,5) )

#%% TEst
#the below example includes plotting and using the metadata of the retrieved data
#import statements
import os
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist

#uncomment the line below, copy data locally and change this path to increase performance
#dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
file_nc_map = os.path.join(dir_testinput,'output','westerscheldt01_map.nc')
file_nc_his = os.path.join(dir_testinput,'output','westerscheldt01_his.nc')

#plot net/grid
ugrid_all = get_netdata(file_nc=file_nc_map)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_aspect('equal')

### test
# from https://gis.stackexchange.com/questions/369026/extract-shapely-polygon-from-matplotlib-hexbin-polycollection

hex_polys = pc.get_paths()[0].vertices
hex_array = []
for xs,ys in pc.get_offsets():
    hex_x = np.add(hex_polys[:,0],  xs)
    hex_y = np.add(hex_polys[:,1],  ys)
    hex_array.append(Polygon(np.vstack([hex_x, hex_y]).T))

hex_grid = gpd.GeoDataFrame({'geometry':hex_array})
gpd.GeoDataFrame({'geometry':hex_array}).plot(edgecolor='black')


data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='s1', timestep=3)#, multipart=False)
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_wl[0,:], ax=None, linewidth=0.5, cmap="jet")
