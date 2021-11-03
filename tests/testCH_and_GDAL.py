# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:12:08 2021

@author: sbe002
Using testConcaveHull-shp

1) gdalgrid to create raster from the bathymetery point
However, a csv from point should be created first along with the .vrt

2) gdalwarp to cut the raster based on the shp (concave)

3) build script to detect max min of no data and delete the tile
"""
#%% This script is from testCreateMatrix
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
from nodataraster import nodataraster

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

matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T

matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])

mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T # create matrix x,y for hull

#% creating concave hull
hull = CH.concaveHull(mat_hull, 9) # utk case ini k=9
# CH.plotPath(matrix, hull)

# create polygon
poly = Polygon(hull)
# plt.plot(*poly.exterior.xy)

import fiona
from fiona.crs import from_epsg
# output location
out_poly_path = 'tests/outputTests/'
out_poly_name = 'CH_'
# Define a polygon feature geometry with one attribute
# project_crs = from_epsg(28992)
schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int', 'area':'float',
                   'length':'float'},
}

# Write a new Shapefile
# source: https://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
with fiona.open((str(out_poly_path)+str(out_poly_name)+'.shp'), 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write({
        'geometry': mapping(poly),
        'properties': {'id': 123, 'area':poly.area,
                       'length':poly.length},
    })

## TODO
# Assign projection after writing shapefile
# Apparently GPKG driver in fiona does not provide polygon
# keep using SHP

#%% Create raster from bathimetri .csv with gdal_grid
# Create .csv file
concave_path = 'tests/outputTests/'
concave_name = 'CH_bathy'
EPSG_coord = 32749

np.savetxt(str(concave_path)+str(concave_name)+'.csv', matrix, fmt="%f", delimiter=",", header='Lon,Lat,Alt',comments='')

# Create .vrt file based on the created .csv
lines = ['<OGRVRTDataSource>', '<OGRVRTLayer name='+'"'+str(concave_name)+'"'+'>',
         '<SrcDataSource>'+str(concave_name)+".csv"+'</SrcDataSource>',
         '<GeometryType>wkbPoint</GeometryType>',
         '<LayerSRS>EPSG:'+str(EPSG_coord)+'</LayerSRS>',
         '<GeometryField separator=" " encoding="PointFromColumns" x="Lon" y="Lat" z="Alt"/>',
         '</OGRVRTLayer>',
         '</OGRVRTDataSource>']
with open(str(concave_path)+str(concave_name)+'.vrt', 'w') as f:
    f.write('\n'.join(lines))
    
# Call and run gdal_grid to create raster file from bathymetry
abs_path = os.path.join('D:/Git/d3d_meso/')
os.chdir(abs_path+str(concave_path))
vrt_in = str(concave_name)+'.vrt'
csv_in = str(concave_name)
ras_out = str(concave_name)+'.tif'
x_min = np.min(matrix[:,0])
x_max = np.max(matrix[:,0])
y_min = np.min(matrix[:,1])
y_max = np.max(matrix[:,1])

x_res = 10
y_res = 10

no_data_val = -999

# create a raster with linear interpolation     
command = 'gdal_grid -zfield "Alt" -a linear:radius=0.0:nodata={no_data_val} -ot Float32 \
            -txe {x_min} {x_max} -tye {y_min} {y_max} \
            -tr {x_res} {y_res} -l {csv_in} {vrt_in} {ras_out} \
            --config GDAL_NUM_THREADS ALL_CPUS'

os.system(command.format(no_data_val=no_data_val, x_min=x_min, x_max=x_max, 
                         y_min=y_min, y_max=y_max, x_res=x_res, y_res=y_res, 
                         csv_in=csv_in, vrt_in=vrt_in, ras_out=ras_out))

# call and run gdalwarp to clip the raster and add no data value
cut_cl = str(out_poly_name)
cut_call = str(out_poly_name)+'.shp'
ras_clip = str(concave_name)+'_clipped'+'.tif'

command_warp = 'gdalwarp -overwrite -of GTiff -cutline {cut_call} -cl {cut_cl} \
                -crop_to_cutline -dstnodata {no_data_val} {ras_out} {ras_clip}'
os.system(command_warp.format(cut_call=cut_call, cut_cl=cut_cl, no_data_val=no_data_val, 
                              ras_out=ras_out, ras_clip=ras_clip))

#%% Integrate raster tiling with IF to delete nodata tiled raster
out_path = "D:/Git/d3d_meso/tests/outputTests/"
output_filename = "tile_"

tile_size_x = 20000 # it is now in meter or pixel if we use script (the last one)
tile_size_y = 20000 # which reads info in pixel

## Set and arrange projection
# Convert the projection with gdal warp
# prep_conv_out = str(out_path) + str(os.path.splitext(input_filename)[0]) + "_UTM" + ".tif"   

ds = gdal.Open(ras_clip)
gt = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
# stats = band.GetStatistics(True,True) # results = min, max, mean, StdDev
# stats[1]-stats[0] = max-min
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
        command = 'gdal_translate -of GTIFF -srcwin {i}, {j}, {xsize_tile}, \
            {ysize_tile} {prep_conv_out} {prep_out}'
        os.system(command.format(i=i, j=j, xsize_tile=xsize_tile, ysize_tile=ysize_tile, 
                                 prep_conv_out=ras_clip, prep_out=prep_out))
        # below is to filter and delete the grid with all nodata value
        nodataraster(prep_out)

