# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 18:02:41 2021

@author: sbe002

### NOTE: big No NO
### Regularized has strange regular pattern
Task for today:
    1. Use scatter_to_regulargrid from dfm_tools to create a regular grid
    2. Instead of saving the file as a netCDF convert that as raster (modify testCH and GDAL)
    3. Take nc2regularGrid as idea and testSplitD3D as a basis of number of grid
    4. Convert the regularised grid as raster via gdal_grid an clip

use the same scripts as in testCH_and_GDAL until line 113
"""
from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid
import os
import numpy as np
from netCDF4 import Dataset
# method in scipy.interpolate.griddata : 'linear','nearest','cubic'




file_nc_map = os.path.join(r'D:\IDM\Tutorial Data - Delft3D FM Suite 2022.01\Tutorial_D-Flow_FM\New_UGrid\Project1.dsproj_data\FlowFM\output\FlowFM_map.nc')
size_per_cell = 200 #m
data_frommap_x = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_x') # this mesh2d_face_x is in the new UGrid
data_frommap_y = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_y')
x_min = np.min(data_frommap_x)
x_max = np.max(data_frommap_x)
y_min = np.min(data_frommap_y)
y_max = np.max(data_frommap_y)

nx = int((x_max-x_min)/ size_per_cell)
ny = int((y_max-y_min)/ size_per_cell)
treg = 50
# tt_val = np.arange(10, 200, 25) # 'all' # can be set as all
lay_val = 'None'
dir_output = 'output2'
key_values = ['time', 'mesh2d_s1', 'mesh2d_face_x', 'mesh2d_face_y'] # pick the one that only required

if not os.path.exists(dir_output):
        os.makedirs(dir_output)
file_nc = file_nc_map
input_nc = Dataset(file_nc, 'r', format='NetCDF4')

time_old = input_nc.variables['time'][:]
if treg != 'all':
    time_old = np.take(time_old, treg)

vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
df = vars_pd
# key_values = ['mesh2d_tem1','time', 'mesh2d_s1', 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY', 'mesh2d_face_x', 'mesh2d_face_y']
df = df.loc[df['nc_varkeys'].isin(key_values)]
data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x') # this mesh2d_face_x is in the new UGrid
data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y') # this one as well

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
        lon = x_grid[0, :]
        lat = y_grid[:, 0]

# import matplotlib.pyplot as plt
field_array2d = np.reshape(field_array, (field_array.shape[0]*field_array.shape[1], field_array.shape[2]))
# plt.contourf(x_grid,y_grid,field_array2d)

y_grid2 = y_grid.T

positions = np.vstack([x_grid.ravel(), y_grid2.ravel(), field_array2d.ravel()]).T

concave_path = 'tests/outputTests/'
concave_name = 'CH_bathy_regular'
np.savetxt(str(concave_path)+str(concave_name)+'.csv', positions, fmt="%f", delimiter=",", header='Lon,Lat,Alt',comments='')

# Call and run gdal_grid to create raster file from bathymetry

abs_path = os.path.join('D:/Git/d3d_meso/')
os.chdir(abs_path+str(concave_path))
vrt_in = str(concave_name)+'.vrt'
csv_in = str(concave_name)
ras_out = str(concave_name)+'.tif'
x_min = np.min(data_frommap_x)
x_max = np.max(data_frommap_x)
y_min = np.min(data_frommap_y)
y_max = np.max(data_frommap_y)

x_res = 10
y_res = 10
     
command = 'gdal_grid -zfield "Alt" -a invdist:power=2.0:smoothing=1.0:nodata=-999 \
            -txe {x_min} {x_max} -tye {y_min} {y_max} \
            -tr {x_res} {y_res} -l {csv_in} {vrt_in} {ras_out} \
            --config GDAL_NUM_THREADS ALL_CPUS'
# command = 'gdal_grid -zfield "Alt" -txe {x_min} {x_max} -tye {y_min} {y_max} \
#             -tr {x_res} {y_res} -l {csv_in} {vrt_in} {ras_out} \
#             --config GDAL_NUM_THREADS ALL_CPUS'
os.system(command.format(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, x_res=x_res, y_res=y_res, 
                         csv_in=csv_in, vrt_in=vrt_in, ras_out=ras_out))

# call and run gdalwarp to clip the raster and add no data value
cut_cl = str(out_poly_name)
cut_call = str(out_poly_name)+'.shp'
ras_clip = str(concave_name)+'_clipped'+'.tif'
no_data_val = -999

command_warp = 'gdalwarp -overwrite -of GTiff -cutline {cut_call} -cl {cut_cl} \
                -crop_to_cutline -dstnodata {no_data_val} {ras_out} {ras_clip}'
os.system(command_warp.format(cut_call=cut_call, cut_cl=cut_cl, no_data_val=no_data_val, 
                              ras_out=ras_out, ras_clip=ras_clip))