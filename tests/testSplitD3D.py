# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:37:17 2021

@author: sbe002
"""

"""
Notes: A succesful splitting test has been attempted.
This scipt has resulted to a consistent splitted raster.
However, the tile size is based on the pixel not metric
"""


## This is the script to test the splitting process
# example: https://howtoinqgis.wordpress.com/2016/12/17/how-to-split-a-raster-in-several-tiles-using-qgis-or-python-gdal/

# or watch this if we want to split the raster into equal pieces
# https://www.youtube.com/watch?v=H5uQ85VXttg

# Import Packages
import os
# import os.path
from osgeo import gdal
# import gdal

# define the file settings
in_path = "D:/Data Topografi/DEMNAS/"
input_filename = "DEMNAS_1016-54_v1.0.tif"
prep_input = str(os.path.join(in_path,input_filename)) #Because path has space trick is made in for loop
prep_input_s = '"' + prep_input + '"'

out_path = "D:/Git/d3d_meso/tests/outputTests/"
output_filename = "tile_"

tile_size_x = 200 # it is now in meter
tile_size_y = 200 

## Set and arrange projection
# Convert the projection with gdal warp
prep_conv_out = str(out_path) + str(os.path.splitext(input_filename)[0]) + "_UTM" + ".tif"   
source_coord = 'EPSG:4326' # if in geographic
target_coord = 'EPSG:32648' #UTM WGS84 48N
com_conv = 'gdalwarp -s_srs {source_coord} -t_srs {target_coord} -r near -of GTiff {prep_input_s} {prep_conv_out}'
os.system(com_conv.format(source_coord=source_coord, target_coord=target_coord, prep_input_s=prep_input_s, prep_conv_out=prep_conv_out))

# =============================================================================
# # assign the projection if no projection is detected
# target_assign_coord = 'EPSG:4326' # if in geographic
# prep_conv_outt = str(out_path) + str(os.path.splitext(input_filename)[0]) + "_assgCoord" + ".tif"  
# # com_as_coord = 'gdal_edit -a_srs {target_assign_coord} {prep_input_s}'
# # com_as_coord = 'gdal_transform -s_srs {target_assign_coord} -t_srs {target_assign_coord}'
# com_as_coord = 'gdalwarp -s_srs {target_assign_coord} -r near -of GTiff {prep_input_s} {prep_conv_outt}'
# os.system(com_as_coord.format(target_assign_coord=target_assign_coord, prep_input_s=prep_input_s, prep_conv_outt=prep_conv_outt))
# =============================================================================


# set the parameters
ds = gdal.Open(in_path + input_filename)
band = ds.GetRasterBand(1)
xsize = band.XSize
ysize = band.YSize

## Trial
ds = gdal.Open(prep_conv_out)
gt = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
xsize = band.XSize
ysize = band.YSize
# get coordinates of upper left corner
xmin = gt[0]
ymax = gt[3]
res = gt[1]

# determine total length of raster
xlen = res * ds.RasterXSize
ylen = res * ds.RasterYSize

# size of a single tile
xsize_tile = int(tile_size_x/res) #num of pixels in tile_size_x
ysize_tile = int(tile_size_y/res) ##num of pixels in tile_size_y

# ----------------------------------

# Tile the raster domain as the prefered tile size in meter
for i in range(0, xsize, xsize_tile):
    for j in range(0, ysize, ysize_tile):
        prep_out = str(out_path) + str(output_filename) + str(i) + "_" + str(j) + ".tif"        
        command = 'gdal_translate -of GTIFF -srcwin {i}, {j}, {xsize_tile}, {ysize_tile} {prep_conv_out} {prep_out}'
        os.system(command.format(i=i, j=j, xsize_tile=xsize_tile, ysize_tile=ysize_tile, prep_conv_out=prep_conv_out, prep_out=prep_out))


# =============================================================================
# # old script with non direct command with os.system
# for i in range(0, xsize, tile_size_x):
#     for j in range(0, ysize, tile_size_y):
#         com_string = "gdal_translate -of GTIFF -srcwin " + str(i) + ", " + \
#         str(j) + ", " + str(tile_size_x) + ", " + str(tile_size_y) + " " + \
#         '"' + prep_input + '"' + " " + str(out_path) + \
#         str(output_filename) + str(i) + "_" + str(j) + ".tif"
#         os.system(com_string)
# =============================================================================

# =============================================================================
# # Tile the raster domain as the prefered tile size in pixel
# for i in range(0, xsize, tile_size_x):
#     for j in range(0, ysize, tile_size_y):
#         prep_out = str(out_path) + str(output_filename) + str(i) + "_" + str(j) + ".tif"        
#         command = 'gdal_translate -of GTIFF -srcwin {i}, {j}, {tile_size_x}, {tile_size_y} {prep_input_s} {prep_out}'
#         os.system(command.format(i=i, j=j, tile_size_x=tile_size_x, tile_size_y=tile_size_y, prep_input_s=prep_input_s, prep_out=prep_out))
# 
# =============================================================================

