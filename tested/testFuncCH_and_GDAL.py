# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:12:08 2021

@author: sbe002

This is the preparation script that is able to:
    1. Read the nc file and create ConcaveHull to derive the outer boundary 
        of the domain 
    2. The ConcaveHull shape as Poly is converted to SHP
    3. The bathi file or other data file (velocity, water level) can be converted
        as raster and clipped based on number 2
    4. The raster can be tiled to particular size
"""
#%% Import packages and set the file
import os
import numpy as np
import numpy.ma as ma
from dfm_tools.get_nc import get_netdata

# import matplotlib.pyplot as plt

import sys
print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles

#uncomment the line below, copy data locally and change this path to increase performance
#dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
# dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
file_nc_map = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_map.nc')
file_nc_his = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_his.nc')
file_nc_ori = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\input\westernscheldt04_net.nc')

#%% Import the nc file, create hull, and build poly
nc_in = file_nc_ori
k_hull = 9
polya = d3dConcaveHull(nc_in,k_hull)

# plt.plot(*polya.exterior.xy)    

#%% From Poly create SHP
projection = 32749
dir_out = os.path.join(r'D:/Git/d3d_meso/tests/outputTests/')
out_poly_name = 'CH_'

d3dPolySHP(polya, dir_out, out_poly_name, projection)

#%% Create raster from bathimetri .csv with gdal_grid
# Create .csv file
concave_path = dir_out
concave_name = 'CH_bathy_'
EPSG_coord = 32749
x_res = 10
y_res = 10
no_data_val = -999
shp_clip = out_poly_name # this assume that the shp file similar with one
                        # build from d3dPolySHP
                        # It will directly refer to the shapefile in the same folder
                        # with 'CH_.shp'
affix = '_clipped'
# create matrix, in this case we can create a custom matrix from data reading
# for instance, selection of data points of WoO
# don't forget to adjust the shapefile from the allowable WoO

ugrid_all = get_netdata(file_nc=file_nc_ori)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])

d3dCSV2ClippedRaster(concave_path, concave_name, EPSG_coord, matrix, x_res, y_res, no_data_val, shp_clip, affix)

#%% Tile the raster and filter + delete nodata tiled raster
out_path = dir_out
output_filename = "tile_"

ras_clip = os.path.join(out_path+concave_name+affix+'.tif')

tile_size_x = 20000 # it is in meter 
tile_size_y = 20000 # which reads info in pixel

d3dRaster2Tiles(out_path, output_filename, ras_clip, tile_size_x, tile_size_y,CreateSHP=True)

import gc
gc.collect() # to clear memory of variables in python after doing del(variables)
