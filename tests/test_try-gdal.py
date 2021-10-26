# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 22:54:51 2021

@author: 
    
    https://help.marine.copernicus.eu/en/articles/5029956-how-to-convert-netcdf-to-geotiff
"""

# packages 
import xarray as xr 
import rioxarray as rio 
import os

# The xarray package enables to load and open the NetCDF file to convert :
yourfile = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01\input\westerscheldt_net.nc')
nc_file = xr.open_dataset(yourfile)

# The next step is to extract the variable of your choice to convert into raster file.
# I want here the bottom temperature bottomT :
bT = nc_file['mesh2d_node_z']

# It is suggested to provide spatial axis x and y for the GeoTIFF 
# and check for Coordinate Reference System (CRS). No output is required.
bT = bT.rio.set_spatial_dims('mesh2d_node_x', 'mesh2d_node_y')
bT.rio.crs

# Copernicus Marine products have as standard CRS the WGS 84 (EPSG 4326) 
# - except for ARCTIC products - however, if necessary, 
# it is possible to define the projection system :
#(Optional)
bT.rio.set_crs("epsg:4326")

# Finally, save the GeoTIFF file : 
bT.rio.to_raster(r"medsea_bottomT_raster.tiff")

#%% 2nd Use Case - GDAL 

# If you are not comfortable with Python or you simply need the .tiff file 
# for your study, a second option is to use the Geospatial Data Abstraction Library GDAL.

# Once installed, you may execute the following command, directly in your terminal prompt : 
gdal_translate NETCDF:"INput_file.nc":variable_name OUTput_filename.tiff

# You will obtain a single-band .tiff containing the variable of interest.
# Here an example of temperature : 
gdal_translate NETCDF:"medsea_tem_20210320_21.nc":thetao med_thetao.tiff

#%% Coba dari notebook
import netCDF4
import os

yourfile = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01\input\westerscheldt_net.nc')
grid =netCDF4.Dataset(yourfile)
# see what variables are in there
grid.variables

# load coordinates
# the [:,:] are required as grid.variables['var'] is only a pointer
lon = grid.variables['mesh2d_node_x'][:,:]
lat = grid.variables['lat'][:,:]

#%% 
from dfm_tools.ugrid import UGrid
    
ugrid = UGrid.fromfile(yourfile)
ugrid.verts.shape[0]

nc_file = xr.open_dataset(yourfile)
nc