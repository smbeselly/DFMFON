# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 19:05:21 2022

@author: sbe002

This is the script for testing the coupling mode

TODO belum nambah script untuk tambahan sediment karena biomass (daun jatuh, dll)
"""
#%% Import the necessary packages, set the file path, and input files

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
# from osgeo import ogr

import geopandas as gpd
import pandas as pd
from xlrd import open_workbook
from xlutils.copy import copy
import shutil
from pathlib import Path
import zipfile
import send2trash

## Set the paths for dll-files and input-files for DFM
PROJ_HOME = os.path.join(r'D:\Git\d3d_meso')
D3D_HOME = os.path.join(r'C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64')
MFON_HOME = os.path.join(PROJ_HOME,'Model-Execute','MesoFON')
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM','FunnelMorphMF30') # model funnel with morpho
MFON_JAR = os.path.join(MFON_HOME, 'complete_model.jar')
gdal_path = os.path.join(r'D:\Program_Files\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils')

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
if not os.path.exists(MFON_Exchange):
    os.makedirs(MFON_Exchange)
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
if not os.path.exists(MFON_Env):
    os.makedirs(MFON_Env)
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
if not os.path.exists(MFON_Trees):
    os.makedirs(MFON_Trees)

dir_out = os.path.join(MFON_Exchange, 'Initialization')
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
## settings

EPSG_Project = 32749 # EPSG code for WGS84/ UTM Zone 49S (Porong case study)
netcdf_domain = os.path.join(D3D_workdir, 'dflowfm', 'Delta_Schematized_funnel_net.nc')
x_res = 10
y_res = 10
no_data_val = -999
tile_size_x = 200 # it is in meter 
tile_size_y = 200 # which reads info in pixel
species_name = 'Avicennia_marina'
LLWL = -1.2

#%% Read the domain and prepare the 'world' for MesoFON
##############
#%% Import the nc file, create hull, and build poly
nc_in = netcdf_domain
LLWL = LLWL - 0.5 # to make correction for concavehull process

ugrid_all = get_netdata(file_nc=nc_in)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# For instance the LLWL value is - 1.5 + 0.5m
matrix_llwl = np.where(matrix[:,2] < LLWL , -999, matrix[:,2])
matrix = np.column_stack([matrix[:,:2],matrix_llwl])
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])
# create matrix x,y for hull
mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T 

# Observe the data first in QGIS and see what kind of concave hull that can be
# created whether an automatic Concave Hull or manual delineation

k_hull = 3 # the lowest number is 3,, you can do the trial and error depend on the model
import ConcaveHull as CH
polya = CH.concaveHull(mat_hull, k_hull)

# plt.plot(*polya.exterior.xy)    

# If concave hull method does not functioning well save the mathull and
# manually process the delineation
np.savetxt(str(dir_out)+ '\\mat_hull'+'.csv', mat_hull, delimiter=",", header='Lon,Lat',comments='')

# df = pd.read_csv(r'Model-Exchange/MesoFON-Env/Tiling/mat_hull_edit.csv')
# mat_hull_df = df.to_numpy()
# polya = CH.concaveHull(mat_hull, k_hull)

#%% From Poly create SHP
## It is needed if we use the ConcaveHull method,
# Otherwise, skip this part and use the manual delineation boundary
# and name the file as CH_.shp
projection = EPSG_Project
dir_out = dir_out
out_poly_name = 'CH_'

d3dPolySHP(polya, dir_out, out_poly_name, projection)

#%% Create raster from bathimetri .csv with gdal_grid
# Create .csv file
concave_path = dir_out
concave_name = 'CH_bathy_'
EPSG_coord = EPSG_Project
x_res = x_res
y_res = y_res
no_data_val = no_data_val
shp_clip = out_poly_name # this assume that the shp file similar with one
                        # build from d3dPolySHP
                        # It will directly refer to the shapefile in the same folder
                        # with 'CH_.shp'
affix = '_clipped'
# create matrix, in this case we can create a custom matrix from data reading
# for instance, selection of data points of WoO
# don't forget to adjust the shapefile from the allowable WoO

ugrid_all = get_netdata(file_nc=nc_in)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])

d3dCSV2ClippedRaster(concave_path, concave_name, EPSG_coord, matrix, x_res, y_res, no_data_val, shp_clip, affix)

#return to home
os.chdir(PROJ_HOME)
#%% Tile the raster and filter + delete nodata tiled raster
out_path = dir_out
output_filename = "tile_"

ras_clip = os.path.join(out_path+str('\\')+concave_name+affix+'.tif')

tile_x = tile_size_x # it is in meter 
tile_y = tile_size_y # which reads info in pixel

d3dRaster2Tiles(out_path, output_filename, ras_clip, tile_x, tile_y,CreateSHP=True)

import gc
gc.collect() # to clear memory of variables in python after doing del(variables)

#%% This part is to tile the mangrove trees file

shp_source = os.path.join(MFON_Trees, 'Master-Trees') #source of the complete trees location
if not os.path.exists(shp_source):
    os.makedirs(shp_source)
shp_source = os.path.join(shp_source, 'MangroveAgeMerged.shp') # name of the trees shp

folder_loc = dir_out # location of the master tiles
file_tile = os.path.join(dir_out,'tile_*.shp' )
save_tiled_trees = os.path.join(MFON_Trees,'Initiate-Trees') # location to save the tiled shp
if not os.path.exists(save_tiled_trees):
    os.makedirs(save_tiled_trees)
    
for filepath in glob.iglob(file_tile):
    # print(filepath)
    save_name = filepath[len(folder_loc):-4]+'_trees'+'.shp'
    # save_loc = os.path.join(folder_loc,save_name)
    save_loc = save_tiled_trees+save_name
    command_ogr = 'ogr2ogr -clipsrc {filepath} {save_loc} {shp_source} -f "ESRI Shapefile"'
    os.system(command_ogr.format(filepath=filepath, save_loc=save_loc, 
                                 shp_source=shp_source))

#%% This part is to write xls file
file_tile_trees = os.path.join(save_tiled_trees,'tile_*_trees.shp')

for filepath in glob.iglob(file_tile_trees): # looping for all with trees affix
    # print(filepath)
    tile_0_read = gpd.read_file(filepath)
    posX = tile_0_read['coord_x']
    posY = tile_0_read['coord_y']
    height_m = tile_0_read['height'] #change for independent files or standardized the shp
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
    xls_name = filepath[len(save_tiled_trees):-4]+'_input'+'.xls'
    xls_loc = save_tiled_trees + xls_name
    
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
    s.write(1,1, species_name)
    s.write(2,0, '/def')
    s.write(3,0, 'values')
        # position_last = start the data + 1 (since it contains header) 
    #  + length of data
    position_last = startrowis + 1 + len(posX)
    s.write(position_last,0, '/values')
    wb.save(xls_loc)

#%% Copy and Place, and Prepare the Tiled xls to the designated folder

file_tile_xls = os.path.join(save_tiled_trees,'tile_*_input.xls')

# baca nama dengan tiles_*_trees.xls
# read ada berapa banyak tile_*_trees
# buat folder masing2 tile di baawh Model-Execute
# Sebelumnya cek apakah ada folder atau tidak, hapus folder lama dan isinya
# Buat folder baru dengan sesuai nama tiles.
# Beri nama MEsoFON_tile_0_0.. dst
# Copy Jar dan extract 
# Read batch_params.xml dan unrolledParamFile dan rubah konten sesuai dengan 
# file kita
# File Raster Salinity di gisSalinityPath
# File Raster Survival (WoO) di gisSurvivalPath

xlsCounter = len(glob.glob(file_tile_xls))

# import pathlib
# print(pathlib.Path(filepath).suffix) # only to get the extension
# recursive delete
# try:
#     shutil.rmtree(os.path.join(MFON_Trees,'tile_*'))
# except OSError as e:
#     print("Error: %s : %s" % (dir_path, e.strerror))

# # Copy the master jar file into each folder
# import shutil

# original = r'C:\Users\Ron\Desktop\Test_1\products.csv'
# target = r'C:\Users\Ron\Desktop\Test_2\products.csv'

# shutil.copyfile(original, target)

### delete old folders, create new tile folders, and extract jar
# delete (recycle bin) the existing tile_* folder
# create the empty folders, copy the jar file, and extract

for filepate in glob.iglob(os.path.join(MFON_HOME,'tile_*')):
    # shutil.rmtree(filepath) # delete recursive
    # print(filepath)
    try:
        send2trash.send2trash(filepate)
    except OSError as e:
        print("Error: %s : %s" % (filepate, e.strerror))

# Create new folders as the existing xls files and extract MesoFON jar file

for filepath in glob.iglob(file_tile_xls): # looping for all with trees affix
    # print(filepath)
   xls_folder = Path(filepath).stem
   os.makedirs(os.path.join(MFON_HOME,xls_folder))
   with zipfile.ZipFile(MFON_JAR, 'r') as zip_ref:
       zip_ref.extractall(os.path.join(MFON_HOME,xls_folder))
 

### Create ideal rasters for initialization files
# 1. delete all content in folder initiate-rasters
# 2. create sal and surv rasters, make copies and give 0 and 1 suffix
# 
    
save_tiled_env = os.path.join(MFON_Env,'Initiate-Rasters') # location to save the tiled shp
try:
    send2trash.send2trash(save_tiled_env)
except OSError as e:
    print("Error: %s : %s" % (save_tiled_env, e.strerror))
  
if not os.path.exists(save_tiled_env):
    os.makedirs(save_tiled_env)


# because gdal_calc is python script it is not straightforward
# source: https://gis.stackexchange.com/questions/268041/inserting-python-variables-into-gdal-calc-expression
gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
# set ideal parameter for initialization (the first 3 month)
calc_surv = '"0*(A>-999)+1"'
calc_sal = '"0*(A>-999)+35"'

for filepatd in glob.iglob(os.path.join(dir_out,'tile_*.tif')):
    # print(filepatd)
    ori_raster = filepatd
    tif_name = Path(filepatd).stem
    surv_raster = os.path.join(save_tiled_env, tif_name+'_surv_.tif')
    sal_raster = os.path.join(save_tiled_env, tif_name+'_sal_.tif')
    # define the syntax
    command_surv = 'python {gdal_calc_path} -A {ori_raster} --outfile={surv_raster} --calc={calc_surv}'
    command_sal = 'python {gdal_calc_path} -A {ori_raster} --outfile={sal_raster} --calc={calc_sal}'
    # calculate
    os.system(command_surv.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
                                  surv_raster=surv_raster, calc_surv=calc_surv))
    os.system(command_sal.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
                                 sal_raster=sal_raster, calc_sal=calc_sal))
    # copy and rename the new raster for surv
    raster_surv = surv_raster
    raster_surv_name = Path(raster_surv).stem
    target_ras_surv0 = os.path.join(save_tiled_env,raster_surv_name+'0.tif')
    target_ras_surv1 = os.path.join(save_tiled_env,raster_surv_name+'1.tif')
    # do the copy-paste
    shutil.copyfile(raster_surv, target_ras_surv0)
    shutil.copyfile(raster_surv, target_ras_surv1)
    # copy and rename the new raster for sal
    raster_sal = sal_raster
    raster_sal_name = Path(raster_sal).stem
    target_ras_sal0 = os.path.join(save_tiled_env,raster_sal_name+'0.tif')
    target_ras_sal1 = os.path.join(save_tiled_env,raster_sal_name+'1.tif')
    # do the copy-paste
    shutil.copyfile(raster_sal, target_ras_sal0)
    shutil.copyfile(raster_sal, target_ras_sal1)
    

### TODO: Find a way how to replace the content in Unrolled_Param and batch_params.xml
# TODO: Ganti rumus yg di line 191, pakai rumus d130










