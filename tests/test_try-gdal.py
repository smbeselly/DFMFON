# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 22:54:51 2021

@author: 
    
    https://help.marine.copernicus.eu/en/articles/5029956-how-to-convert-netcdf-to-geotiff
"""

# packages 
import xarray as xr 
import rioxarray as rio 

# The xarray package enables to load and open the NetCDF file to convert :
nc_file = xr.open_dataset('./mydirectory/medsea_tem_20200320_21.nc')

# The next step is to extract the variable of your choice to convert into raster file.
# I want here the bottom temperature bottomT :
bT = nc_file['bottomT']

# It is suggested to provide spatial axis x and y for the GeoTIFF 
# and check for Coordinate Reference System (CRS). No output is required.
bT = bT.rio.set_spatial_dims('lon', 'lat')
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