# -*- coding: utf-8 -*-
"""
Copyright notice
-------------------------
This library is developed as part of the PhD research of
Sebrian Mirdeklis Beselly Putra conducted at IHE Delft Institute 
for Water Education and Delft University of Technology

The author  of this library is:
    Sebrian Beselly
    s.besellyputra@un-ihe.org
    s.m.beselly@tudelft.nl
    sebrian@ub.ac.id
    
    IHE Delft Institute for Water Education,
    PO Box 3015, 2601DA Delft
    the Netherlands
    
This library is free software: you can redistribute it and/or modify 
it under the terms of the GPL-3.0 license
    
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GPL-3.0 license for more details.

Publication related to this library
Beselly, S.M., U. Grueters, M. van Der Wegen, J. Reyns, J. Dijkstra, and D. Roelvink. “Modelling Mangrove-Mudflat Dynamics with a Coupled Individual-Based-Hydro-Morphodynamic Model.” Environmental Modelling & Software, August 28, 2023, 105814. https://doi.org/10.1016/j.envsoft.2023.105814.


@author: sbe002
source: https://github.com/openearth/dflowfm_regularize_output/blob/master/nc2regularGrid_listComprehension.py
This is a modified version of the OpenEarht by FineWilms

Modifications:
    1. add the 'if function' to filter the function for 3D or 2D model, indicated by the layers information
    2. add dir_output option, so that user can freely decide where to store the regularized version

Description:
    This function creates a regularized grid of the new D3D-FM UGrid
    Output of the function is a regular grid with the specified x and y
    and value of the water level (mesh2d_s1) and water velocity in form of
    velocity in x and y (mesh2d_ucx and mesh2d_ucy)
    
Parameters:
    fileNC : string, is the path to the 0'th partition of the DFlowFM map file output
    xpts : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
    ypts : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
    tms : numpy array or 'all', an array of times you want to do the interpolation for
    lrs : numpy array, 'all', or integer. The number of layers you wish to include. The script detect if there are layers or not. 
    dir_output: specify where to save the output
    
    It accepts the value of 
    tms = 'all', or integer of time (25), or as np.array
        # tms = np.arange(0,4,2)
    lrs = 'all' for 3D, 'None' for 2D, or as np.array
        # lrs = np.arange(0,3) 
    dir_output : can be relative to the cwd, example: 'output'
        or abs path with os
    

Example:
    fpin = os.path.join(Simulation, 'Output', 'Flow_map.nc')
    nx = 400
    ny = 500
    treg = np.arange(10, 200, 25) #or 'all' or 25
    lreg = 'None' #because it is 2D or 'all' or np.arange(0,3) in 3D
    dir_output = 'output_reg'
    
    regularGrid_to_netcdf(fp_in, nx, ny, treg, lreg, dir_output):

"""

from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid
import os
import numpy as np
from netCDF4 import Dataset
# import time as tm

def regularGrid_to_netcdf(fp_in, nx, ny, treg, lreg, dir_output):
    # dir_output = os.path.abspath(os.path.join(os.path.dirname(__file__),'..', 'output'))
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    file_nc = fp_in
    input_nc = Dataset(file_nc, 'r', format='NetCDF4')
    time_old = input_nc.variables['time'][:]
    if treg != 'all':
        time_old = np.take(time_old, treg)

    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)


    df = vars_pd

    key_values = ['mesh2d_tem1','time', 'mesh2d_s1', 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY', 'mesh2d_face_x', 'mesh2d_face_y']

    df = df.loc[df['nc_varkeys'].isin(key_values)]

    """
    ####################################################################################################################
    #   Regularise all files with 3 dimensions (time, nFaces, layers). 
    #   This will be equal to four dimensions in the regular grid format since nFaces is the x- and y- dimension.
    ####################################################################################################################
    """
    # if lreg == 'None': #make sure 2D model is loaded without layers
    #     df2 = df.loc[df['ndims'] == 2] # I don't know why but the model has ndims = 2 in original script ndims==3
    # else:
    #     df2 = df.loc[df['ndims'] == 3]
      
    data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x') # this mesh2d_face_x is in the new UGrid
    data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y') # this mesh2d_face_y is in the new UGrid
    # time = get_ncmodeldata(file_nc=file_nc, varname='time', timestep=treg)
    # outname = '%s_regular.nc' % os.path.split(fileNC)[1][0:-3]
    outname = '%s_regular.nc' % os.path.split(file_nc)[1][0:-3]
    file_nc_reg = os.path.join(dir_output, outname)
    root_grp = Dataset(file_nc_reg, 'w', format='NETCDF4')
    # root_grp = Dataset(file_nc_reg, 'w', file_format='NETCDF4')
    root_grp.description = 'Regularized the map output'
    first_read = True
    i = 0
    
    if lreg != 'None':
       
        df2 = df.loc[df['ndims'] == 3]  
       
        for index, row in df2.iterrows(): # I changed from df to df2
    
            if row['dimensions'][1] == 'mesh2d_nEdges':
                continue
            
            data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row['nc_varkeys'], timestep=treg, layer=lreg)    
            data_frommap_var = data_frommap_var.filled(np.nan)
            field_array = np.empty((data_frommap_var.shape[0], ny, nx, data_frommap_var.shape[-1]))
            tms = data_frommap_var.shape[0]
            lrs = data_frommap_var.shape[-1]
            trange = range(0, tms)
            lrange = range(0, lrs)
    
            A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                                 values=data_frommap_var[t, :, l].flatten(), method='linear') for t in
                          trange for l in lrange])
            x_grid = A[0][0]
            y_grid = A[0][1]
            A = A[:, 2, :, :]
            A = np.moveaxis(A, [0], [2])
            subs = np.split(A, tms, axis=2)
    
            field_array[:, :, :, 0:lrs] = [subs[tn] for tn in trange]
            field_array = np.ma.masked_invalid(field_array)
            print('done with variable %s' % row['nc_varkeys'])
    
            if first_read:
                unout = 'seconds since 2015-01-01 00:00:00'
                lon = x_grid[0, :]
                lat = y_grid[:, 0]
                # create dimensions
                root_grp.createDimension('time', None)
                root_grp.createDimension('lon', lon.shape[0])
                root_grp.createDimension('lat', lat.shape[0])
                root_grp.createDimension('layer', lrs)
                lonvar = root_grp.createVariable('lon', 'float32', 'lon')
                lonvar.setncattr('axis', 'X')
                lonvar.setncattr('reference', 'geographical coordinates, WGS84 projection')
                lonvar.setncattr('units', 'degrees_east')
                lonvar.setncattr('_CoordinateAxisType', 'Lon')
                lonvar.setncattr('long_name', 'longitude')
                lonvar.setncattr('valid_max', '180')
                lonvar.setncattr('valid_min', '-180')
                lonvar[:] = lon
    
                latvar = root_grp.createVariable('lat', 'float32', 'lat')
                latvar.setncattr('axis', 'Y')
                latvar.setncattr('reference', 'geographical coordinates, WGS84 projection')
                latvar.setncattr('units', 'degrees_north')
                latvar.setncattr('_CoordinateAxisType', 'Lat')
                latvar.setncattr('long_name', 'latitude')
                latvar.setncattr('valid_max', '90')
                latvar.setncattr('valid_min', '-90')
                latvar[:] = lat
    
                layervar = root_grp.createVariable('layer', 'float32', 'layer')
                layervar.setncattr('axis', 'Z')
                layervar.setncattr('reference', 'geographical coordinates, WGS84 projection')
                layervar.setncattr('units', 'm')
                layervar.setncattr('_CoordinateZisPositive', 'down')
                layervar.setncattr('_CoordinateAxisType', 'Height')
                layervar.setncattr('long_name', 'Depth')
    
                layervar[:] = range(0, lrs)
    
                timevar = root_grp.createVariable('time', 'float32', 'time')
                timevar.setncattr('units', unout)
                timevar.setncattr('calendar', 'standard')
                timevar.setncattr('long_name', 'time')
                timevar.setncattr('_CoordinateAxisType', 'Time')
                
                timevar[:] = time_old
    
            fieldName = row['nc_varkeys']
            fieldvar = root_grp.createVariable(fieldName, 'float32', ('time', 'lat', 'lon', 'layer'), fill_value=-999)
            key = fieldName
            for ncattr in input_nc.variables[key].ncattrs():
                if ncattr != "_FillValue":
                    root_grp.variables[fieldName].setncattr(ncattr, input_nc.variables[key].getncattr(ncattr))
    
    
            fieldvar[:] = field_array
            first_read = False
            i += 1


    """
    ####################################################################################################################
    #   Regularise all files with 2 dimensions (time, nFaces, layers).
    #   This will be equal to 3 dimensions in the regular grid format since nFaces is the x- and y- dimension.
    ####################################################################################################################
    """
    # else:
    if lreg == 'None':
        print('STARTING 2D')
    df2 = df.loc[df['ndims'] == 2]
    
    excludeList = ['edge', 'face', 'x', 'y']
    for index, row in df2.iterrows():
        test = any(n in str(row['nc_varkeys']) for n in excludeList)
        if not test:
            if row['dimensions'][1] == 'mesh2d_nEdges':
                continue
            ntimes = row['shape'][0]
            data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row['nc_varkeys'], timestep=treg)
            data_frommap_var = data_frommap_var.filled(np.nan)
            field_array = np.empty((data_frommap_var.shape[0], ny, nx))
            trange = range(0, data_frommap_var.shape[0])
            tms = data_frommap_var.shape[0]
            A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny, 
                                                 values=data_frommap_var[t, :].flatten(), method='linear') for t in trange])
    
            x_grid = A[0][0]
            y_grid = A[0][1]
            A = A[:, 2, :, :] # select the parameters
            field_array[:, :, :] = A
            field_array = np.ma.masked_invalid(field_array)
            # field_array = np.ma.masked_invalid(A)
            
            """"set netCDF parameters"""
            
            unout = str(df.loc[df['nc_varkeys'] == 'time']['units'].values)[2:-9]
                  
            lon = x_grid[0, :]
            lat = y_grid[:, 0]
            # create dimensions
            root_grp.createDimension('time', None)
            root_grp.createDimension('lon', lon.shape[0])
            root_grp.createDimension('lat', lat.shape[0])
            lonvar = root_grp.createVariable('lon', 'float32', 'lon')
            lonvar.setncattr('axis', 'X')
            # lonvar.setncattr('reference', 'geographical coordinates, WGS84 projection')
            # lonvar.setncattr('units', 'degrees_east')
            lonvar.setncattr('_CoordinateAxisType', 'Lon')
            lonvar.setncattr('long_name', 'longitude')
            # lonvar.setncattr('valid_max', '180')
            # lonvar.setncattr('valid_min', '-180')
            lonvar[:] = lon
    
            latvar = root_grp.createVariable('lat', 'float32', 'lat')
            latvar.setncattr('axis', 'Y')
            # latvar.setncattr('reference', 'geographical coordinates, WGS84 projection')
            # latvar.setncattr('units', 'degrees_north')
            latvar.setncattr('_CoordinateAxisType', 'Lat')
            latvar.setncattr('long_name', 'latitude')
            # latvar.setncattr('valid_max', '90')
            # latvar.setncattr('valid_min', '-90')
            latvar[:] = lat
      
            timevar = root_grp.createVariable('time', 'float32', 'time')
            timevar.setncattr('units', unout)
            timevar.setncattr('calendar', 'standard')
            timevar.setncattr('long_name', 'time')
            timevar.setncattr('_CoordinateAxisType', 'Time')
      
            timevar[:] = time_old
    
            """write data to new netcdf"""
            fieldName = row['nc_varkeys']
            fieldvar = root_grp.createVariable(fieldName, 'float32', ('time', 'lat', 'lon'), fill_value=-999)
            key = fieldName
            for ncattr in input_nc.variables[key].ncattrs():
                if ncattr != "_FillValue":
                    root_grp.variables[fieldName].setncattr(ncattr, input_nc.variables[key].getncattr(ncattr))
            fieldvar[:] = field_array
    root_grp.close()
