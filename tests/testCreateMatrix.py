# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 10:57:29 2021

@author: sbe002
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
    
# %%Translate np array to raster
### Still not yet get the nice result
import rasterio as rio    
prep_out = os.path.join(r'D:/Git/d3d_meso/tests/outputTests')
data_matrix = ma.compressed(matrix)
with rio.open('outname.tif', 'w', driver='GTiff') as dst:
    dst.write(matrix, 1)

def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


def main(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):
    reversed_arr = array[::-1] # reverse array so the tif looks like the array
    array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,reversed_arr) # convert array to raster

rasterOrigin = (np.min(matrix[:,0]),np.min(matrix[:,1]))
pixelWidth = 10
pixelHeight = 10
newRasterfn = 'test.tif'

array = np.array([[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])

main(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,matrix)

#%% rasterio
# source: https://www.youtube.com/watch?v=bZMgVe6s33s
x = np.linspace(-4.0, 4.0, 240)
y = np.linspace(-3.0, 3.0, 180)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-2 * np.log(2) * ((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 1 ** 2)
Z2 = np.exp(-3 * np.log(2) * ((X + 0.5) ** 2 + (Y + 0.5) ** 2) / 2.5 ** 2)
Z = 10.0 * (Z2 - Z1)
plt.imshow(Z)
plt.contour(Z)

from rasterio.transform import Affine
res = (x[-1] - x[0]) / 240.0
transform = Affine.translation(x[0] - res / 2, y[0] - res / 2) * Affine.scale(res, res)

with rio.open(
    'new.tif',
    'w',
    driver='GTiff',
    height=Z.shape[0],
    width=Z.shape[1],
    count=1,
    dtype=Z.dtype,
    crs='+proj=latlong',
    # transform=transform,
) as dst:
    dst.write(Z, 1)
    
#%%
# Create 3D Matrix
from scipy.interpolate import griddata
data_x = ma.compressed(ugrid_all.mesh2d_node_x)
data_y = ma.compressed(ugrid_all.mesh2d_node_y)
data_z = ma.compressed(ugrid_all.mesh2d_node_z)

points = np.block([[data_x],[data_y]]).T

#Assign the raster resolution
rasterRes = 100

# get the dimension
xDim = np.max(points[:,0])-np.min(points[:,0])
yDim = np.max(points[:,1])-np.min(points[:,1])
nCols = int(xDim/rasterRes) #num of pixels in tile_size_x
nRows = int(yDim/rasterRes) ##num of pixels in tile_size_y

#%% Test GDAL
import numpy as np
from osgeo import gdal, gdal_array, osr
import matplotlib.pylab as plt

array = np.array(( (0.1, 0.2, 0.3, 0.4),
                   (0.2, 0.3, 0.4, 0.5),
                   (0.3, 0.4, 0.5, 0.6),
                   (0.4, 0.5, 0.6, 0.7),
                   (0.5, 0.6, 0.7, 0.8) ))
# My image array      
lat = np.array(( (10.0, 10.0, 10.0, 10.0),
                 ( 9.5,  9.5,  9.5,  9.5),
                 ( 9.0,  9.0,  9.0,  9.0),
                 ( 8.5,  8.5,  8.5,  8.5),
                 ( 8.0,  8.0,  8.0,  8.0) ))
lon = np.array(( (20.0, 20.5, 21.0, 21.5),
                 (20.0, 20.5, 21.0, 21.5),
                 (20.0, 20.5, 21.0, 21.5),
                 (20.0, 20.5, 21.0, 21.5),
                 (20.0, 20.5, 21.0, 21.5) ))
# For each pixel I know it's latitude and longitude.
# As you'll see below you only really need the coordinates of
# one corner, and the resolution of the file.

# Source: https://gis.stackexchange.com/questions/37238/writing-numpy-array-to-raster-file
xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
nrows,ncols = np.shape(array)
xres = (xmax-xmin)/float(ncols)
yres = (ymax-ymin)/float(nrows)
geotransform=(xmin,xres,0,ymax,0, -yres)   
# That's (top left x, w-e pixel resolution, rotation (0 if North is up), 
#         top left y, rotation (0 if North is up), n-s pixel resolution)
# I don't know why rotation is in twice???

output_raster = gdal.GetDriverByName('GTiff').Create('myraster.tif',ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
srs = osr.SpatialReference()                 # Establish its coordinate encoding
srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.
                                             # Anyone know how to specify the 
                                             # IAU2000:49900 Mars encoding?
output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
                                                   # to the file
output_raster.GetRasterBand(1).WriteArray(array)   # Writes my array to the raster

output_raster.FlushCache()

# x = np.linspace(0, 2, 6)
# y = np.linspace(0,1,6)
# z = np.linspace(0,3,6)
# X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
# plt.imshow(Z)
# plt.contour(Z)


# x_ = np.linspace(0., 1., 10)
# y_ = np.linspace(1., 2., 20)
# z_ = np.linspace(3., 4., 30)

# x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

# a = np.block([[lon[:,0]],[lat[:,0]],[array[:,0]]]).T
# np.savetxt("foo.csv", a, delimiter=",")

# gridOpt = gdal.GridOptions(format='GTiff', zfield=array[:,0])
# output = gdal.Grid('outcome.tif', [lon[:,0],lat[:,0]])
# option: https://gdal.org/python/index.html
# source : https://www.youtube.com/watch?v=FSnJ2VXNV3c
# source https://www.qgistutorials.com/en/docs/interpolating_point_data.html
a = np.block([[lon[:,0]],[lat[:,0]]]).T
idw = gdal.Grid("invdist.tif", a, zfield=array[:,0],
                outputBounds=[ulx,uly,lrx,lry],
                width=xsize, height=ysize)

"""
 TODO 
 mungkin caranya adalah interpolasi dengan IDW (misal), kemudian di-clip 
 sesuai dengan master
 Kemudian di-buat logika bahwa semua nilai di bawah atau di atas nilai
 tertentu pada raster adalah equal dengan atau high salinity
"""