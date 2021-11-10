# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 18:43:08 2021

@author: sbe002
This script is able to:
Define the approximated area for mangrove to grow (LLWL+0.5m) --> DONE
Tile the trees location (.shp) similar with the tiles in raster and convert
    as (.xls) input files

"""

# =============================================================================
# To support the 2.5 assignment, a new function capabilities in d3d_prep_raster 
# has been created.
# It will create tiles of shp for the raster by define CreateSHP as True
# =============================================================================

## To support 2.5 a script by adapting the d3dRaster2Tiles by only taking the shp
# creation can be used#

import os
import numpy as np
import numpy.ma as ma

import sys
print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles
from dfm_tools.get_nc import get_netdata
import matplotlib.pyplot as plt
plt.close('all')

import glob
from osgeo import ogr

import geopandas as gpd
import pandas as pd
from xlrd import open_workbook
from xlutils.copy import copy

#%% This is to create a ConcaveHull for the filtered depth of the allowable
#   mangrove area

#uncomment the line below, copy data locally and change this path to increase performance
# file_nc_map = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_map.nc')
# file_nc_his = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_his.nc')
file_nc_ori = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\input\westernscheldt04_net.nc')
LLWL = -1.5
LLWL = LLWL - 0.5 # to make correction for concavehull process

ugrid_all = get_netdata(file_nc=file_nc_ori)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# For instance the LLWL value is - 1.5 + 0.5m
matrix_llwl = np.where(matrix[:,2] < -2 , -999, matrix[:,2])
matrix = np.column_stack([matrix[:,:2],matrix_llwl])
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])
# create matrix x,y for hull
mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T 

# Observe the data first in QGIS and see what kind of concave hull that can be
# created whether an automatic Concave Hull or manual delineation

# =============================================================================
# # Use this for automatic delineation
# import ConcaveHull as CH
# hull = CH.concaveHull(mat_hull, 3)
# # create polygon
# from shapely.geometry import Polygon
# poly = Polygon(hull)
# plt.plot(*poly.exterior.xy)
# =============================================================================

#%% This part is to tile the mangrove trees file

shp_source = os.path.join(r'D:/Git/d3d_meso/tests/inputTests/Temp_For_Trial.shp')
folder_loc = os.path.join(r'D:/Git/d3d_meso/tests/outputTests')
file_tile = os.path.join(folder_loc,'tile_*.shp' )

for filepath in glob.iglob(file_tile):
    # print(filepath)
    save_name = filepath[len(folder_loc):-4]+'_trees'+'.shp'
    # save_loc = os.path.join(folder_loc,save_name)
    save_loc = folder_loc+save_name
    command_ogr = 'ogr2ogr -clipsrc {filepath} {save_loc} {shp_source} -f "ESRI Shapefile"'
    os.system(command_ogr.format(filepath=filepath, save_loc=save_loc, 
                                 shp_source=shp_source))
    # print(file_tile)

#%% This part is to write xls file
# def totuple(a):
#     try:
#         return tuple(totuple(i) for i in a)
#     except TypeError:
#         return a

# ### Read Shapefile
# tile_0_0 = os.path.join(r'D:/Git/d3d_meso/tests/inputTests/North ttops edited.shp')

# tile_0_read = gpd.read_file(tile_0_0)
# coordx = tile_0_read['coord x']
# coordy = tile_0_read['coord y']
# coordz = tile_0_read['CHM_North']

# coba = np.block([[coordx],[coordy],[coordz]]).T
# coba2 = totuple(coba)

# ex: https://xlsxwriter.readthedocs.io/tutorial01.html

###### Trial with real dataset

file_tile_trees = os.path.join(folder_loc,'tile_*_trees.shp' )

for filepath in glob.iglob(file_tile_trees): # looping for all with trees affix
    # print(filepath)
    tile_0_read = gpd.read_file(filepath)
    posX = tile_0_read['coord x']
    posY = tile_0_read['coord y']
    height_m = tile_0_read['CHM_North']
    id_id = np.arange(len(posX))
    speciesName = np.ones(len(posX))*1 #if only one species is recorded, however it is better to place this in shapefile
    types_species = id_id+1 # this variable starts from 1 to N+1
    shiftedBelowPos = np.ones(len(posX))*1
    age = np.ones(len(posX))*1 # should be derived from shapefile #TO DO will be added later
    rbh_m = (0.0015*height_m**2)+(0.0015*height_m) #dummy it uses equation in test_trial
    
    # Create Panda Dataframe with Header
    df = pd.DataFrame({'id': id_id,
                       'speciesName' : speciesName,
                       'type' : types_species,
                       'posX' : posX,
                       'posY' : posY,
                       'shiftedPosX' : posX,
                       'shiftedPosY' : posY,
                       'shiftedBelowPosX' : speciesName,
                       'shiftedBelowPosY' : speciesName,
                       'age' : age,
                       'rbh_m' : rbh_m,
                       'newRbh_m' : rbh_m,
                       'height_m' : height_m,
                       'newHeight_m' : height_m})
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    xls_name = filepath[len(folder_loc):-4]+'_input'+'.xls'
    xls_loc = folder_loc + xls_name
    
    # Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
    # the default header and skip one row to allow us to insert a user defined
    # header.
    startrowis = 4
    
    with pd.ExcelWriter(xls_loc) as writer:
        df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=startrowis, header=True)
    
    rb = open_workbook(xls_loc)
    wb = copy(rb)
    # Insert the value in worksheet
    s = wb.get_sheet(0)
    s.write(0,0,'def')
    s.write(1,0, 1)
    s.write(1,1, 'Rhizophora_apiculata')
    s.write(2,0, '/def')
    s.write(3,0, 'values')
        # position_last = start the data + 1 (since it contains header) 
    #  + length of data
    position_last = startrowis + 1 + len(posX)
    s.write(position_last,0, '/values')
    wb.save(xls_loc)

