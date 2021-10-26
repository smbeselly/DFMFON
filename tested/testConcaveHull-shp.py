# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:14:02 2021

@author: sbe002

This script can create a concave hull and write shapefile based on the 
concave
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
from shapely.geometry import mapping, Polygon
import geopandas as gpd
from osgeo import gdal, ogr, os, osr

import sys
print(sys.path)
sys.path.append('FnD3D')
import ConcaveHull as CH

#uncomment the line below, copy data locally and change this path to increase performance
#dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
file_nc_map = os.path.join(dir_testinput,'output','westerscheldt01_map.nc')
file_nc_his = os.path.join(dir_testinput,'output','westerscheldt01_his.nc')
file_nc_ori = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01\input\westerscheldt_net.nc')

#get all ugrid net data
ugrid_all = get_netdata(file_nc=file_nc_ori)#,multipart=False)
# data_x = ugrid_all.mesh2d_node_x
# data_y = ugrid_all.mesh2d_node_y
# data_z = ugrid_all.mesh2d_node_z
# plt.scatter(data_x,data_y)

# unmask the data
# data_xx = ma.compressed(data_x)

# matrix = np.block([[data_x],[data_y],[data_z]]).T

matrix = np.block([[ugrid_all.mesh2d_node_x],[ugrid_all.mesh2d_node_y],[ugrid_all.mesh2d_node_z]]).T

matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])

mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T # create matrix x,y for hull
#%% creating concave hull
hull = CH.concaveHull(mat_hull, 5) # utk case ini k=9
CH.plotPath(matrix, hull)

# create polygon
poly = Polygon(hull)
plt.plot(*poly.exterior.xy)

import fiona
from fiona.crs import from_epsg

# Define a polygon feature geometry with one attribute
project_crs = from_epsg(28992)
schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int', 'area':'float',
                   'length':'float'},
}

# Write a new Shapefile
# source: https://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
with fiona.open('my_shp2.shp', 'w', 'ESRI Shapefile', project_crs, schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write({
        'geometry': mapping(poly),
        'properties': {'id': 123, 'area':poly.area,
                       'length':poly.length},
    })