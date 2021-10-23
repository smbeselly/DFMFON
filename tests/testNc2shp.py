# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 21:09:55 2021

@author: sbe002
This is a trial script to convert netcdf to shp
"""

import netCDF4
import numpy as np
from osgeo import gdal,osr,ogr
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import xarray as xr
import os


# read in file path for shapefile
fp_shp = os.path.join("D:/Git/d3d_meso/tests/cb_2018_us_ua10_500k/cb_2018_us_ua10_500k.shp")
# read in netcdf file path
dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
ncs = os.path.join(dir_testinput,'output','westerscheldt01_map.nc')


# Read in NETCDF as a pandas dataframe
# Xarray provides a simple method of opening netCDF files, and converting them to pandas dataframes
ds = xr.open_dataset(ncs)
edgar = ds.to_dataframe()

# the index in the df is a Pandas.MultiIndex. To reset it, use df.reset_index()
edgar = edgar.reset_index()

# Read shapefile using gpd.read_file()
shp = gpd.read_file(fp_shp)

# read the netcdf data file
# nc = netCDF4.Dataset(ncs,'r')

# quick check for shpfile plotting
shp.plot(figsize=(12, 8));

# filter out shapefile for SPECIFIC city/region

# how to filter rows in DataFrame that contains string
# extract NYC from shapefile dataframe
nyc_shp = shp[shp['NAME10'].str.contains("New York")]

# export shapefile
#nyc_shp.to_file('NYC.shp', driver ='ESRI Shapefile')

# use geopandas points_from_xy() to transform Longitude and Latitude into a list of shapely.Point objects and set it as a geometry while creating the GeoDataFrame
edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.lon, edgar.lat))

print(edgar_gdf.head())

# check CRS coordinates
nyc_shp.crs #shapefile
edgar_gdf.crs #geodataframe netcdf

# set coordinates equal to each other
# PointsGeodataframe.crs = PolygonsGeodataframe.crs
edgar_gdf.crs = nyc_shp.crs

# check coordinates after setting coordinates equal to each other
edgar_gdf.crs #geodataframe netcdf

# Clip points, lines, or polygon geometries to the mask extent.
mask = gpd.clip(edgar_gdf, nyc_shp)

