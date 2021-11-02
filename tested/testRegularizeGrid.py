# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 13:27:48 2021

@author: sbe002

Aplication of nc2regulargrid for processing of ouput data
source: https://github.com/openearth/dflowfm_regularize_output/blob/master/nc2regularGrid_listComprehension.py
"""

import sys
# print(sys.path)
sys.path.append('FnD3D')
import os
import numpy as np

from nc2regularGrid_listComprehension_mod import regularGrid_to_netcdf


# dir_testinput = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01')
# file_nc_map = os.path.join(dir_testinput,'output','westerscheldt01_map.nc')
# file_nc_his = os.path.join(dir_testinput,'output','westerscheldt01_his.nc')
# file_nc_ori = os.path.join(r'D:\Delft3D FM\Tutorial_D-Flow_FM\tutorial09\tutorial09.dsproj_data\westerscheldt01\input\westerscheldt_net.nc')

file_nc_map = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_map.nc')
x_val = 500
y_val = 200
# tt_val = 50
tt_val = np.arange(10, 200, 25) # 'all' # can be set as all
lay_val = 'None'
dir_output = 'output2'

"""fileNC : string, is the path to the 0'th partition of the DFlowFM map file output
    xpts : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
    ypts : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
    tms : numpy array or 'all', an array of times you want to do the interpolation for
    lrs : numpy array, 'all', or integer. The number of layers you wish to include. The script detect if there are layers or not. 
    """

regu = regularGrid_to_netcdf(file_nc_map, x_val, y_val, tt_val, lay_val, dir_output)