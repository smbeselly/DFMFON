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
"""

#%% Input Folders and Files
PROJ_HOME = r'E:\06.RepastDFM\Simulation_Results\Research_Run\d3d_meso_run10_template_rev\d3d_meso_run10_scenario_A'
SYS_APP = r'E:\06.RepastDFM\Simulation_Results\Research_Run\d3d_meso_run10_template_rev\d3d_meso_run10_scenario_A\FnD3D'
D3D_HOME = r'F:\Research_Run\DFM_EXE\2sebrian_20220518' #sal_veg-OK
gdal_loc = r'C:\Users\brian\anaconda3\envs\d3dfm\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\brian\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'
Mangr_SHP = 'Tip_Saplings_geser_2.shp'

D3D_Model = 'd3d_meso_run10_scenario_A'
D3D_Domain = 'Grid_Funnel_1_net.nc'
MFON_Folder = 'MesoFON_20220506_noSurv'

## Check the complete_model.jar file and change this source file
Sal_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_sal_'
Surv_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_surv_'
Excel_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_trees_input.xls'

#%% Import the necessary packages, set the file path, and input files

import os
import numpy as np
import numpy.ma as ma

import sys
print(sys.path)
sys.path.append(SYS_APP) # as this Func will be in the same folder, no longer needed
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles#, D3dNewRaster2Tiles
# from d3d_meso_mangro import csv2ClippedRaster
from dfm_tools.get_nc import get_netdata
from d3d_meso_mangro import calcAgeCoupling0
import matplotlib.pyplot as plt
plt.close('all')

import glob
# from osgeo import ogr
# Supress/hide the warning
np.seterr(invalid='ignore')

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d
from xlrd import open_workbook
from xlutils.copy import copy
import shutil
from pathlib import Path
from osgeo import gdal, gdalconst
import zipfile
import send2trash
import sys
import fileinput
import re

# to log the print to file : https://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("logfile_01_prep.log", "a")
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass    

sys.stdout = Logger()

## Set the paths for dll-files and input-files for DFM
PROJ_HOME = os.path.join(PROJ_HOME)
D3D_HOME = os.path.join(D3D_HOME)
MFON_HOME = os.path.join(PROJ_HOME,'Model-Execute','MesoFON')
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM',D3D_Model) 
MFON_JAR = os.path.join(MFON_HOME, MFON_Folder,'complete_model.jar')
MFON_LocalBatchRunner = os.path.join(MFON_HOME,'local_batch_run.properties')
gdal_path = os.path.join(gdal_loc)
JAVAREP = os.path.join(JAVA_Exe)

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
if not os.path.exists(MFON_Exchange):
    os.makedirs(MFON_Exchange)
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
if not os.path.exists(MFON_Env):
    os.makedirs(MFON_Env)
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
if not os.path.exists(MFON_Trees):
    os.makedirs(MFON_Trees)
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')
if not os.path.exists(MFON_OUT):
    os.makedirs(MFON_OUT)


dir_out = os.path.join(MFON_Exchange, 'Initialization')
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
## settings

EPSG_Project = 32749 # EPSG code for WGS84/ UTM Zone 49S (Porong case study)
netcdf_domain = os.path.join(D3D_workdir, 'dflowfm', D3D_Domain)
x_res = 10
y_res = 10
no_data_val = -999.0
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

### Observe the data first in QGIS and see what kind of concave hull that can be

# df = pd.read_csv(r'Model-Exchange/MesoFON-Env/Tiling/mat_hull_edit.csv')
# mat_hull_df = df.to_numpy()
# polya = CH.concaveHull(mat_hull, k_hull)

#%% From Poly create SHP
## It is needed if we use the ConcaveHull method,
# Otherwise, skip this part and use the manual delineation boundary
# and name the file as CH_.shp
# If you are a QGIS user, use Concave hull (k-nearest neighbor)
# set k=3 as start

projection = EPSG_Project
dir_out = dir_out
out_poly_name = 'CH_'

# =============================================================================
# d3dPolySHP(polya, dir_out, out_poly_name, projection)
# =============================================================================

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
shp_source = os.path.join(shp_source, Mangr_SHP) # name of the trees shp

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
    age = tile_0_read['age']

    if height_m.size != 0:
        height_m = height_m.values  # in metre
        # rbh_m = np.where(height_m < 1.37, f_0(height_m*100)/100/2, f_dbh(height_m*100)/100/2)
        rbh_m = 0.03/2
    else:
        rbh_m = tile_0_read['height']

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

#%% Copy and Place the payload.jar, and Prepare the Tiled xls to the designated folder

file_tile_xls = os.path.join(save_tiled_trees,'tile_*_input.xls')

# Delete the existing the execution folder tile_* to prepare for the new one
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
   shutil.copy(MFON_LocalBatchRunner, os.path.join(MFON_HOME,xls_folder))

 

### Create ideal rasters for initialization files

save_tiled_env = os.path.join(MFON_Env,'Initiate-Rasters') # location to save the tiled shp
try:
    send2trash.send2trash(save_tiled_env)
except OSError as e:
    print("Initiate-Rasters Folder is already deleted before this command: %s : %s" % (save_tiled_env, e.strerror))
  
if not os.path.exists(save_tiled_env):
    os.makedirs(save_tiled_env)

# Use gdal_calc.py to generate the ideal rasters for intializing the MesoFON
# because gdal_calc is python script, therefore it is not straightforward
# source: https://gis.stackexchange.com/questions/268041/inserting-python-variables-into-gdal-calc-expression
gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
# set ideal parameter for initialization (the first 3 month)
calc_surv = '"0*(A>-999)+1"'
calc_sal = '"0*(A>-999)+35"'
val_no_data_surv = 0
val_no_data_sal = 60

def changeNoData(datatif,value):
    maskfile = gdal.Open(datatif, gdalconst.GA_Update)
    maskraster = maskfile.ReadAsArray()
    maskraster = np.where((maskraster > -999), maskraster, value ) # yg lebih dari -999 diganti jadi 0
    maskband = maskfile.GetRasterBand(1)
    maskband.WriteArray( maskraster )
    maskband.FlushCache()

for filepatd in glob.iglob(os.path.join(dir_out,'tile_*.tif')):
    # print(filepatd)
    ori_raster = filepatd
    tif_name = Path(filepatd).stem
    surv_raster = os.path.join(save_tiled_env, tif_name+'_surv_.tif')
    sal_raster = os.path.join(save_tiled_env, tif_name+'_sal_.tif')
    # define the syntax
    command_surv = 'python {gdal_calc_path} -A {ori_raster} --outfile={surv_raster} --calc={calc_surv} --NoDataValue={no_data_val}'
    command_sal = 'python {gdal_calc_path} -A {ori_raster} --outfile={sal_raster} --calc={calc_sal} --NoDataValue={no_data_val}'
    # calculate
    os.system(command_surv.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
                                  surv_raster=surv_raster, calc_surv=calc_surv, no_data_val=no_data_val))
    os.system(command_sal.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
                                 sal_raster=sal_raster, calc_sal=calc_sal, no_data_val=no_data_val))
    # change the nodata value into 0 for surv and 60 for sal
    changeNoData(surv_raster,val_no_data_surv)
    changeNoData(sal_raster,val_no_data_sal)
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
    

### Replace the content in Unrolled_Param, batch_params.xml, and parameters.xml
# function to replace the content of the file
# source = https://www.delftstack.com/howto/python/python-replace-line-in-file/
def replacement(file, previousw, nextw):
   for line in fileinput.input(file, inplace=1):
       line = line.replace(previousw, nextw)
       sys.stdout.write(line)
# The Looping to replace the variables in the params files
for filepatf in glob.iglob(os.path.join(MFON_Exchange,'Initialization','tile_*.tif')):
    Init_Rasters = save_tiled_env
    Init_Trees = save_tiled_trees
    # Create Update Path
    Sal_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_sal_').replace("\\",'/')
    Surv_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_surv_').replace("\\",'/')
    Excel_Update = os.path.join(Init_Trees,Path(filepatf).stem+'_trees_input.xls').replace("\\",'/')
    # Point the Files that need to be updated
    unrolledParam = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','unrolledParamFile.txt')
    batchParam = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','scenario.rs','batch_params.xml')
    param_xml = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','scenario.rs','parameters.xml')
    # Replace the content of the file
    # unrolledParamFile.txt
    replacement(unrolledParam, Surv_Source, Surv_Update)       
    replacement(unrolledParam, Sal_Source, Sal_Update)
    replacement(unrolledParam, Excel_Source, Excel_Update)
    # batch_params.xml
    replacement(batchParam, Surv_Source, Surv_Update)       
    replacement(batchParam, Sal_Source, Sal_Update)
    replacement(batchParam, Excel_Source, Excel_Update)
    # parameters.xml
    replacement(param_xml, Surv_Source, Surv_Update)       
    replacement(param_xml, Sal_Source, Sal_Update)
    replacement(param_xml, Excel_Source, Excel_Update)
    
#%% Use os.system(command) to call Java from Python    

for filepatt in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
    if gpd.read_file(os.path.join(save_tiled_trees,Path(filepatt).stem[:-6]+'.shp')).size > 0:
    # only calculate MesoFON if trees exist 
        # delete the existing instance1 in in the folder to prevent symlink errorp prior running
        try:
            send2trash.send2trash(os.path.join(filepatt,'instance_1'))
        except OSError as e:
            print("Instance 1 is already deleted before this command: %s : %s" % (os.path.join(filepatt,'instance_1'), e.strerror))
        
        # cd to the directory where MesoFon Exec is located
        os.chdir(filepatt)
        print('Run MesoFON model', Path(filepatt).stem)
        command_java = '{JAVAREP} -cp lib/* repast.simphony.batch.LocalDriver local_batch_run.properties'
        os.system(command_java.format(JAVAREP=JAVAREP))
        # back to the project home
        os.chdir(PROJ_HOME)


### retrieving the results and copy to the MesoFON Model-Out
namae=[] #empty list for initializing the namae
for filepatg in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
    nama = []
    for name in glob.iglob(os.path.join(filepatg, 'instance_1','MF_Trees_*.txt')):
        nama.append(name)
    nama = list(filter(lambda x: not re.search('batch_param_map', x), nama)) # exclude batch_param.txt
    MFON_OUT_tile = os.path.join(MFON_OUT,Path(filepatg).stem)
    if not os.path.exists(MFON_OUT_tile):
        os.makedirs(MFON_OUT_tile)
    # select the MFON_Trees only and paste it to the MesoFON_Out
    if nama != []:
        shutil.copyfile(nama[0], os.path.join(MFON_OUT_tile,Path(nama[0]).stem + ' Coupling_0.txt'))
        namae.append(nama[0])

### Compile the results to compile folder
MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')
if not os.path.exists(MFON_OUT_compile):
    os.makedirs(MFON_OUT_compile)

all_df = []    
for nama_a in namae:
    df = pd.read_csv(nama_a)
    all_df.append(df)

Concat_table = pd.concat(all_df)
# drop tick 0 year, only take 0.25
Concat_table = Concat_table[Concat_table.tick > 0]
run_is = 'Coupling_0' # change this with the real name

# use spatial in scipy to match the x,y of the mangroves and the age information.
master_trees = gpd.read_file(shp_source)
age_coupling = calcAgeCoupling0(Concat_table, master_trees)
Concat_table['Age'] = age_coupling #update age after MesoFON run
Concat_table.drop(['tick'], inplace=True, axis=1)
Concat_table = Concat_table.reset_index(drop=True) # reset the index

# Concatenated table is saved as txt file
Concat_table.to_csv(os.path.join(MFON_OUT_compile,run_is+'.txt'), sep=',', index=False, header=True)

### retrieving canopy out to the MesoFON Model-Out
namac=[] #empty list for initializing the namae
for filepatg in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
    nama = []
    for name in glob.iglob(os.path.join(filepatg, 'instance_1','MF_Canopy_*.txt')):
        nama.append(name)
    nama = list(filter(lambda x: not re.search('batch_param_map', x), nama)) # exclude batch_param.txt
    MFON_OUT_tile = os.path.join(MFON_OUT,Path(filepatg).stem)
    if not os.path.exists(MFON_OUT_tile):
        os.makedirs(MFON_OUT_tile)
    # select the MFON_Trees only and paste it to the MesoFON_Out
    if nama != []:
        shutil.copyfile(nama[0], os.path.join(MFON_OUT_tile,Path(nama[0]).stem + ' Coupling_0.txt'))
        namac.append(nama[0])

### Compile the results to compile folder
MFON_OUT_compile_can = os.path.join(MFON_OUT,'Compile_Canopy')
if not os.path.exists(MFON_OUT_compile_can):
    os.makedirs(MFON_OUT_compile_can)

all_dfc = []    
for nama_a in namac:
    df = pd.read_csv(nama_a)
    all_dfc.append(df)

Concat_tablec = pd.concat(all_dfc)
Concat_tablec = Concat_table.reset_index(drop=True) # reset the index
# Concatenated table is saved as txt file
Concat_tablec.to_csv(os.path.join(MFON_OUT_compile_can,run_is+'.txt'), sep=',', index=False, header=True)

