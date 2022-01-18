# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 13:50:58 2022

@author: sbe002

This the scipt to do the running after the preparation and
initialization of the MesoFON and DFM model.
"""

# =============================================================================
# 1. Read the compiled text file
# 2. Do the isolate the trees and calculate the drag coefficient 
# 3. set_var the cdveg with new drag coefficient calculated by the algorithm 
# 4. Do it in a loop
# =============================================================================
#%% Reset all variable
# Don't forget to %reset %clear

#%% Import the packages and set the file
import numpy as np
import numpy.ma as ma
import bmi.wrapper
# import ctypes
# import cmocean.cm
# import matplotlib.colors
# import matplotlib.pyplot as plt
# import matplotlib.tri as tri
import os
# import datetime
import pandas as pd
# from scipy import integrate
import faulthandler
faulthandler.enable()
import sys
# print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed

# from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc import get_ncmodeldata
# from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes, calcDragCoeff 
from d3d_meso_mangro import calcWOO, calcAgeCoupling0, createPointSHP, createXLSfromSHP  
from d3d_meso_mangro import createRaster4MesoFON, modifyParamMesoFON, calcDragInLoop

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

import gc
import geopandas as gpd
from pathlib import Path
import glob
# from osgeo import gdal, gdalconst
import send2trash
import re
import shutil


## Set the paths for dll-files and input-files for DFM
PROJ_HOME = os.path.join(r'D:\Git\d3d_meso')
D3D_HOME = os.path.join(r'C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64')
MFON_HOME = os.path.join(PROJ_HOME,'Model-Execute','MesoFON')
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM','FunnelMorphMF30_Adjusted') # model funnel with morpho
MFON_JAR = os.path.join(MFON_HOME, 'complete_model.jar')
MFON_LocalBatchRunner = os.path.join(MFON_HOME,'local_batch_run.properties')
gdal_path = os.path.join(r'D:\Program_Files\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils')
JAVAREP = os.path.join(r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe')

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')

dir_out = os.path.join(MFON_Exchange, 'Initialization')

figsavefolder= os.path.join(PROJ_HOME,'Model-Out','Figures')
if not os.path.exists(figsavefolder):
    os.makedirs(figsavefolder)
#%% Settings
EPSG_Project = 32749 # EPSG code for WGS84/ UTM Zone 49S (Porong case study)
coupling_period = 30 #days, actually 90 days
MorFac = 30 
woo_inun = 3 # inundation free period (days)
species_name = 'Avicennia_marina'

x_res = 10
y_res = 10
no_data_val = -999.0
tile_size_x = 200 # it is in meter 
tile_size_y = 200 # which reads info in pixel
species_name = 'Avicennia_marina'
shp_clip = os.path.join(dir_out, 'CH_shp')

concave_name = 'CH_bathy_'
affix = '_clipped'

## Check the complete_model.jar file and change this source file
# for later to be updated with the new params file as generated during the preprocessing
# Sal_Source = r'F:\Temp\MesoFONbatch_JDK11\Data_Trees_JDK11\Raster_Dummy_UTM_'
# Surv_Source = r'F:\Temp\MesoFONbatch_JDK11\Data_Trees_JDK11\Raster_Dummy_UTM_Surv_'
# Excel_Source = r'F:\Temp\MesoFONbatch_JDK11\Data_Trees_JDK11\test_trial.xls'
Sal_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_sal_'
Surv_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_surv_'
Excel_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_trees_input.xls'


#%% Initiate the BMI
# search and locate the files
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm','bin')
dflowfm_engine = os.path.join(dflowfm_path, 'dflowfm.dll')
config_file = os.path.join(D3D_workdir, 'FunnelMorphMF30_Adjusted.xml') # funnel morpho

mdu_file = os.path.join(D3D_workdir, 'dflowfm', 'FlowFM.mdu')
grid_file = os.path.join(D3D_workdir, 'dflowfm', 'Grid_Funnel_1_net.nc')


### Access the information from mesh
mesh_face_x = get_ncmodeldata(file_nc=grid_file, varname='mesh2d_face_x')
mesh_face_x = ma.compressed(mesh_face_x)
mesh_face_y = get_ncmodeldata(file_nc=grid_file, varname='mesh2d_face_y')
mesh_face_y = ma.compressed(mesh_face_y)
mesh_face_nodes = get_ncmodeldata(file_nc=grid_file, varname='mesh2d_face_nodes')

### Initialize BMI Model
## Add corrects locations to environment variable PATH for DFM
os.environ['PATH'] = os.path.join(D3D_HOME, 'share', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dflowfm', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dimr', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dwaves', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'esmf', 'scripts') \
+ ";" + os.path.join(D3D_HOME, 'swan', 'scripts')

## Define DFM wrapper
# os.chdir(dflowfm_path)
model_dfm = bmi.wrapper.BMIWrapper(engine=dflowfm_engine, configfile=mdu_file)

## Define and initialise DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()

#%% Access the variable and calculate the cell number
# Access the variables from BMI.get_var
xk = model_dfm.get_var('xk') #Net node x coordinate {"shape": ["numk"]}
yk = model_dfm.get_var('yk')
#xz, yz, related with bl, s1, 
xz = model_dfm.get_var('xz') #waterlevel point / cell centre, x-coordinate (m) {"location": "face", "shape": ["ndx"]}
yz = model_dfm.get_var('yz') #y coordinate
#xzw, yzw related with cell number
xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
yzw = model_dfm.get_var('yzw') #y coordinate
# get the cell area

#### calculate the cell number as in the position of xzw and yzw or ndxi
xyzw_cell_number = create_xyzwCellNumber(xzw,yzw,mesh_face_x,mesh_face_y)
xyzw_nodes = create_xyzwNodes(mesh_face_nodes,xyzw_cell_number)

#### Read Master Trees
master_trees = gpd.read_file(os.path.join(MFON_Trees, 'Master-Trees', 'MangroveAgeMerged.shp'))

#%% Read the compiled tree from the MesoFON Initialization Run and calculate drag_coefficient

MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')

read_data = pd.read_csv(os.path.join(MFON_OUT_compile,'Coupling_0.txt'))
# use spatial in scipy to match the x,y of the mangroves and the age information.
# age_coupling0 = calcAgeCoupling0(read_data, master_trees)
age_coupling = calcAgeCoupling0(read_data, master_trees)

# For loop for all of the cell number

drag_coeff = [] #should explore whether as list or as array list

for row in range(len(xyzw_cell_number)):
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    # find the position based on the cell number
    position = xyzw_cell_number[row,2].astype(int)
    
    nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
    nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
    # Find the min max of each x,y coordinate
    # create the list of x_min-x_max and y_min-y_max
    x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
    y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
    
    # subsetting pandas 
    read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                 (read_data['GeoRefPosX'] <= x_range[1])]
    read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                        (read_data_subset['GeoRefPosY'] <= y_range[1])]
    # TODO check cell_area and water_depth
    cell_area = ba[row] # check this later 
    water_depth = hs[row]# adjust
    trees_data = read_data_subset
    # calculate the drag coefficient (currently only for Avicennia marina)
    cd_veg = calcDragCoeff(x_range, y_range, cell_area, water_depth, trees_data)
    # append the calculated cd_veg to the drag_coeff list
    drag_coeff.append(cd_veg)


# update the variable with new value
model_dfm.set_var('Cdvegsp',drag_coeff)

#%% Loop the Coupling

coupling_period = coupling_period*24*3600 # change from days to second
coupling_period_model = coupling_period/MorFac #real use in model (in seconds)
# how many loops required by Delft3D to finish one coupling period
coupling_time = coupling_period_model/model_dfm.get_time_step() 
# how many coupling is required
coupling_ntime = model_dfm.get_end_time()/coupling_period_model
# to accommodate not integer coupling ntime, take the floor value
# and in the end of looping continue the rest of the simulation if the value is
# not integer
coupling_ntimeUse = np.floor(coupling_ntime)
        
for ntime in range(int(coupling_ntimeUse)):
    # do the calculation for each coupling_ntime
    ### 1. run the DFM all the simulation time within ntime
    water_level = np.empty((len(xyzw_cell_number),0)) 
    #https://www.delftstack.com/howto/numpy/python-numpy-empty-array-append/
    for itime in range(int(coupling_time)):
        model_dimr.update()
        s1 = model_dfm.get_var('s1') 
        # store the maximum water level per time step in column wise
        water_level = np.append(water_level, np.reshape(s1,(len(xyzw_cell_number),1)), axis=1) #append the s1
    
    ### 2. convert the water_level from each time_step to each day
    col_num = (24*3600)/ model_dfm.get_time_step() # equal to how many columns represent 1 day
    # find daily maximum water level
    wl_shape = water_level.shape
    h_wl = np.empty((len(xyzw_cell_number),0))
    for ii in range(int(wl_shape[1]/col_num)):
        bb = np.amax(water_level[:,ii*col_num:ii*col_num+col_num])
        h_wl = np.append(h_wl,bb)
    
    ### 3. Calculate the probability based on the WoO and store the value in each cell
    fromheight = [] # empty list
    Pvalue = [] # empty list
    for ii in range(h_wl.shape[0]):
        fromheightcalc, Pvaluecalc = calcWOO(h_wl[ii,:],woo_inun) # get correlation of elevation and Probability
        # Append the value as list
        fromheight = np.append(fromheightcalc) #append as list
        Pvalue = np.append(Pvaluecalc)
    # 3.1. calculate the probability of WoO based on median value
    med_h_wl = np.median(h_wl, axis=1) # find median value of the h_wl for each cell number
    # 3.2. create the survival probability array based on WoO
    surv_val = np.empty(len(med_h_wl)) #initiate an empty array
    for row in range(med_h_wl.shape):
        surv_val[row] = np.interp(med_h_wl[row],fromheight[row],Pvalue[row])
    
    ### 4. Create the raster from the surv-val
    surv_val_raster = np.column_stack((xzw,yzw,surv_val))
    concave_path = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
    if not os.path.exists(concave_path):
        os.makedirs(concave_path)
    d3dCSV2ClippedRaster(concave_path, concave_name, EPSG_Project, surv_val_raster, x_res, y_res, no_data_val, shp_clip, affix)
    #return to home
    os.chdir(PROJ_HOME)
    
    ### 5. Tile the raster and filter + delete nodata tiled raster
    # this is the master of surv tile with tile_0_0 etc...
    output_filename = "tile_"
    ras_clip = os.path.join(concave_path, concave_name+affix+'.tif')
    d3dRaster2Tiles(concave_path, output_filename, ras_clip, tile_size_x, tile_size_y,CreateSHP=True)
    gc.collect() # to clear memory of variables in python after doing del(variables)
    
    ### 6. Create Mangrove Trees shp and tile
    # Master Trees of Each Coupling
    #TODO age_coupling harus diganti biar generik
    # hasil simulasi kok mangrove height-nya besar sekali,, perlu dicek
    # Create master pt shp with info: coord_x, coord_y, height_m, age, rbh_m
    # master tree name: coupling+str(ntime+1)
    # 6.1. createPointSHP(read_data, age_coupling0, concave_path, EPSG_Project) 
    createPointSHP(read_data, age_coupling, concave_path, EPSG_Project) 
    
    # 6.2. Tile the Master Trees
    shp_source = os.path.join(concave_path, 'coupling'+str(ntime+1)) # location of the source tree and the master tree shp
    folder_loc = dir_out # location of the master tiles
    file_tile = os.path.join(dir_out,'tile_*.shp' )
    save_tiled_trees = os.path.join(MFON_Trees,'coupling'+str(ntime+1)) # location to save the tiled shp
    if not os.path.exists(save_tiled_trees):
        os.makedirs(save_tiled_trees)
        
    for filepath in glob.iglob(file_tile):
        # print(filepath)
        save_name = Path(filepath).stem+'_trees'+'.shp'
        # save_loc = os.path.join(folder_loc,save_name)
        save_loc = save_tiled_trees+save_name
        command_ogr = 'ogr2ogr -clipsrc {filepath} {save_loc} {shp_source} -f "ESRI Shapefile"'
        os.system(command_ogr.format(filepath=filepath, save_loc=save_loc, 
                                     shp_source=shp_source))
    
    ### 7. Create xls file of the tiled trees
    # however, it still uses the Avicennia Marina 
    file_tile_trees = os.path.join(save_tiled_trees,'tile_*_trees.shp')
       
    # 7.1. define d_137 equation for Avicennia : taken from Uwe Grueter's calculation
    a0 = -0.172
    b0 = 49.0713765855412
    a137 = -0.172
    b137 = 48.10139
    
    # 7.2. create the XLS files
    createXLSfromSHP(file_tile_trees, a0, b0, a137, b137, save_tiled_trees, species_name)

    ### 8. Create the env raster and tile raster files
    # location of the tile is in concave_path\tile_0_0.tif
    # 1. create new folder: coupling+str(ntimes+1)
    # 2. create sal rasters, make copies and give 0 and 1 suffix 
    # 3. copy the surv rasters as created in Create the raster from the surv-val
    #    make copies and give 0 and 1 suffix 
    save_tiled_env = os.path.join(MFON_Env, 'coupling'+str(ntime+1))
    if not os.path.exists(save_tiled_env):
        os.makedirs(save_tiled_env)
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    # calc_surv = '"0*(A==-999)"'
    calc_sal = '"0*(A>-999)+35"'
    # val_no_data_surv = 0
    val_no_data_sal = 60
            
    createRaster4MesoFON(concave_path, gdal_calc_path, no_data_val, 
                         save_tiled_env, calc_sal, val_no_data_sal)
    
    ### 9. Replace the content in Unrolled_Param, batch_params.xml, and parameters.xml
    # function to replace the content of the file
    
    Surv_Source, Sal_Source, Excel_Source = modifyParamMesoFON(MFON_HOME, MFON_Exchange, 
                                                               save_tiled_env, save_tiled_trees, 
                                                               Surv_Source, Sal_Source, Excel_Source)
    
    ### 10. Run the MesoFON
    for filepatt in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
        # delete the existing instance1 in in the folder to prevent symlink errorp prior running
        try:
            send2trash.send2trash(os.path.join(filepatt,'instance_1'))
        except OSError as e:
            print("Error: %s : %s" % (os.path.join(filepatt,'instance_1'), e.strerror))
        # cd to the directory where MesoFon Exec is located
        os.chdir(filepatt)
        print('Run MesoFON model', Path(filepatt).stem)
        command_java = '{JAVAREP} -cp lib/* repast.simphony.batch.LocalDriver local_batch_run.properties'
        os.system(command_java.format(JAVAREP=JAVAREP))
        # back to the project home
        os.chdir(PROJ_HOME)


    ### 11. retrieving the results and copy to the MesoFON Model-Out
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
        shutil.copyfile(nama[0], os.path.join(MFON_OUT_tile,Path(nama[0]).name))
        namae.append(nama[0])

    ### 12. Compile the results to compile folder
    MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')
    if not os.path.exists(MFON_OUT_compile):
        os.makedirs(MFON_OUT_compile)

    all_df = []    
    for nama_a in namae:
        df = pd.read_csv(nama_a)
        all_df.append(df)

    Concat_table = pd.concat(all_df)
    # 12.1. drop tick 0 year, only take 0.25
    Concat_table = Concat_table[Concat_table.tick > 0]
    run_is = 'coupling'+str(ntime+1) # change this with the real name
    # 12.2. Concatenated table is saved as txt file
    Concat_table.to_csv(os.path.join(MFON_OUT_compile, run_is+'.txt'), sep=',', index=False, header=True)
    
    ### 13. Read the compile txt file and create the Cd
    # read_data = pd.read_csv(os.path.join(MFON_OUT_compile, run_is+'.txt'))
    read_data = Concat_table
    # 13.1. use spatial in scipy to match the x,y of the mangroves and the age information.
    master_trees = os.path.join(concave_path+str('\\')+Path(concave_path).stem+'.shp')
    age_coupling = calcAgeCoupling0(read_data, master_trees) # prepare the age_coupling for the Master Trees of Each Coupling procedure

    # 14. For loop for all of the cell number

    drag_coeff = calcDragInLoop(xyzw_cell_number, model_dfm, xyzw_nodes, xk, yk, read_data)


    # 15. update the variable with new value
    model_dfm.set_var('Cdvegsp',drag_coeff)
    
### End Loop
#Finalize the running
model_dimr.finalize()   


# to check if the value is not integer
# coupling_ntimeUse == coupling_ntime

# test = r"D:\Git\d3d_meso\Model-Exchange\Initialization\tile_0_20.shp"   


















