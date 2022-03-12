# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 13:50:58 2022

@author: sbe002

This is the script for model include salinity and include seedlings establishments
"""

#%% Reset all variable
# Don't forget to %reset %clear

#%% Input Folders and Files
PROJ_HOME = r'D:\Git\d3d_meso'
SYS_APP = r'D:\Git\d3d_meso/FnD3D'
D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\oss_artifacts_x64_140691'
gdal_loc = r'D:\Program_Files\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'

D3D_Model = 'FunnelMorphMF30_Adjusted_Saline_geser_2'
D3D_Domain = 'Grid_Funnel_1_net.nc'
config_xml = 'FunnelMorphMF30_Adjusted_Saline_geser.xml'
mdu_file = 'FlowFM.mdu'

Mangr_SHP = 'geserMangroveAgeMerged.shp'

#%% Import the packages and set the file
import numpy as np
import numpy.ma as ma
import bmi.wrapper
# import matplotlib.path as mpltPath
import os
import pandas as pd
import faulthandler
faulthandler.enable()
import sys
# print(sys.path)
sys.path.append(SYS_APP) # as this Func will be in the same folder, no longer needed

from dfm_tools.get_nc import get_ncmodeldata
# from dfm_tools.io.polygon import Polygon
from dfm_tools.io.mdu import read_deltares_ini
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes #, calcDragCoeff 
from d3d_meso_mangro import calcWOO, calcAgeCoupling0, createPointSHP #, createXLSfromSHP  
from d3d_meso_mangro import modifyParamMesoFON #, createRaster4MesoFON, calcDragInLoop
from d3d_meso_mangro import csv2ClippedRaster, Sald3dNewRaster2Tiles #, clipSHPcreateXLSfromGPD
from d3d_meso_mangro import SalNew_func_createRaster4MesoFON, newCalcDraginLoop
from d3d_meso_mangro import New_clipSHPcreateXLSfromGPD, SalNew_Sal_func_createRaster4MesoFON
from d3d_mangro_seeds import index_veg, seedling_establishment
from d3d_mangro_seeds import seedling_dispersal, calculate_residual, collect_res
from d3d_mangro_seeds import seedling_prob

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE" #to prevent error in matplotlib

import geopandas as gpd
from pathlib import Path
import glob
import send2trash
import re
import shutil
from scipy.interpolate import interp1d
from dateutil import parser
import datetime
# import copy

## Set the paths for dll-files and input-files for DFM
PROJ_HOME = os.path.join(PROJ_HOME)
D3D_HOME = os.path.join(D3D_HOME)
MFON_HOME = os.path.join(PROJ_HOME,'Model-Execute','MesoFON')
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM',D3D_Model) # model funnel with morpho
# MFON_JAR = os.path.join(MFON_HOME, 'complete_model.jar')
MFON_LocalBatchRunner = os.path.join(MFON_HOME,'local_batch_run.properties')
gdal_path = os.path.join(gdal_loc)
JAVAREP = os.path.join(JAVA_Exe)

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')

dir_out = os.path.join(MFON_Exchange, 'Initialization')

figsavefolder= os.path.join(PROJ_HOME,'Model-Out','Figures')
if not os.path.exists(figsavefolder):
    os.makedirs(figsavefolder)
seedlings_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings')
if not os.path.exists(seedlings_out):
    os.makedirs(seedlings_out)
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
shp_clip = os.path.join(dir_out, 'CH_shp')

concave_name = 'CH_bathy_'
affix = '_clipped'

## Check the complete_model.jar file and change this source file
Sal_Source = 'Initiate-Rasters'
Surv_Source = Sal_Source
Excel_Source = 'Initiate-Trees'


#%% Initiate the BMI
# search and locate the files
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm','bin')
dflowfm_engine = os.path.join(dflowfm_path, 'dflowfm.dll')
config_file = os.path.join(D3D_workdir, config_xml) # funnel morpho

mdu_file = os.path.join(D3D_workdir, 'dflowfm', mdu_file)
grid_file = os.path.join(D3D_workdir, 'dflowfm', D3D_Domain)


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

#### calculate the cell number as in the position of xzw and yzw or ndxi
xyzw_cell_number = create_xyzwCellNumber(xzw,yzw,mesh_face_x,mesh_face_y)
xyzw_nodes = create_xyzwNodes(mesh_face_nodes,xyzw_cell_number)

#### Read Master Trees
master_trees = gpd.read_file(os.path.join(MFON_Trees, 'Master-Trees', Mangr_SHP))

#### Read the polygon pli and add the indices
# pli = Polygon.fromfile(os.path.join(D3D_workdir,'dflowfm','vege.pli'))
# path = mpltPath.Path(pli[0][0])
# ind = path.contains_points(np.transpose((xz,yz))).astype(int)

#%% Read the compiled tree from the MesoFON Initialization Run and calculate drag_coefficient

MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')

read_data = pd.read_csv(os.path.join(MFON_OUT_compile,'Coupling_0.txt'))
read_data['Height_cm'] = read_data['Height_cm']

# use spatial in scipy to match the x,y of the mangroves and the age information.
# age_coupling0 = calcAgeCoupling0(read_data, master_trees)
# age_coupling = calcAgeCoupling0(read_data, master_trees)
age_coupling = read_data['Age']

# For loop for all of the cell number
index_veg_cel = index_veg(model_dfm, xyzw_cell_number, xyzw_nodes, xk, yk, read_data)
drag_coeff = newCalcDraginLoop(xyzw_cell_number, xyzw_nodes, xk, yk, 
                               read_data, model_dfm, index_veg_cel)
# assume that the boundary flow nodes are located at the end of array
addition = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + 0.005  
# drag_coeff = np.append(drag_coeff, addition)*ind
drag_coeff = np.append(drag_coeff, addition)

#%% Get time reference from the model

## get the reference date and starttime of model
getmdu = read_deltares_ini(mdu_file)
refdate = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'RefDate')]['value']
tstart = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'TStart')]['value']

# parse dfm's time in string to datetime var in Python
refdatet = parser.parse(refdate.iloc[0])
tstartt = datetime.timedelta(seconds=float(tstart.iloc[0]))
#reference time is
reftime = refdatet+tstartt

#%% Loop the Coupling
# change from days to second
coupling_period = coupling_period*24*3600 
# time required in model to achieve same coupling period (in seconds)
coupling_period_model = coupling_period/MorFac 
# how many loops required by Delft3D to finish one coupling period or
# how many iterations is needed to achieve the coupling period model
coupling_time = coupling_period_model/model_dfm.get_time_step() 
# how many coupling is required or needed with MesoFON
coupling_ntime = model_dfm.get_end_time()/coupling_period_model # how many coupling is 
# to accommodate not integer coupling ntime, take the floor value
# and in the end of looping continue the rest of the simulation if the value is
# not integer
# if the number is not round use the floor value
coupling_ntimeUse = np.floor(coupling_ntime) 
curyr_check = 0
list_seed2sapl = []
# add_seeds_age = coupling_period / (24*3600)
add_seeds_age = datetime.timedelta(days = coupling_period / (24*3600))
        
for ntime in range(int(coupling_ntimeUse)):
    # do the calculation for each coupling_ntime
    print('Start the coupling',str(ntime+1),'computation')
    ### 1.1. run the DFM all the simulation time within ntime
    # update the variable with new value taken from previous drag calculation
    model_dfm.set_var('Cdvegsp',drag_coeff)
    water_level = np.empty((len(xz),0)) 
    salinity = np.empty((len(xz),0)) 
    #https://www.delftstack.com/howto/numpy/python-numpy-empty-array-append/
    
    # find cells that have vegetation
    index_veg_cel = index_veg(model_dfm, xyzw_cell_number, xyzw_nodes, xk, yk, read_data)
        
    try:
        list_seed2sapl.append(seedling_finalpos)
        seed2sapl = pd.concat(list_seed2sapl, axis=0)
        # add age to the appended pandas
        seed2sapl['Age'] = seed2sapl['Age']+add_seeds_age
        
        list_seed2sapl = []
        list_seed2sapl.append(seed2sapl)
        
        check_as_saplings = seed2sapl[seed2sapl['Age'] >= datetime.timedelta(days = 730)]

        if check_as_saplings.size > 0:
            print('Transfer seedlings as saplings at coupling', str(ntime+1))
            # select seedlings that have been transformed to saplings
            now_as_saplings = check_as_saplings.copy()
            now_as_saplings.loc[:, ['Age']] = (now_as_saplings['Age']/datetime.timedelta(days=365))
            now_as_saplings['dbh_cm'] = 3
            now_as_saplings['Height_cm'] = 100
            
            now_as_saplings.drop(['ParentPosX', 'ParentPosY','row',
                                  'Created Time Stamp', 'CrownSurfaceArea_m2'], 
                                 inplace=True, axis=1)
            
            # update the read data with new saplings as mature mangroves
            read_data = pd.concat([read_data, now_as_saplings], axis=0, ignore_index=True)
            
            # reset list to use the filterd list from current selection
            # filter for less than 730
            filter_seeds2apl = seed2sapl[seed2sapl['Age'] <= datetime.timedelta(days = 730)]
            # list_seed2sapl = []
            seed2sapl = filter_seeds2apl.copy()
            list_seed2sapl = []
            list_seed2sapl.append(seed2sapl)
            
            print('Total number of saplings:', now_as_saplings.shape[0])
            
        else:
            print("The seedlings' age is less than 2 years")
    except:
        print('seedlings production is not yet initiated (1st run)')

    
    ### 1.2. Check for seedling establishment
    # seedling establishment only occur during fruiting season (January)
    # current simulation time is
    cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) #https://stackoverflow.com/questions/775049/how-do-i-convert-seconds-to-hours-minutes-and-seconds
    # current simulation time multiply by MF
    cursecMF = cursec*MorFac
    # current month based on simulation time (with MF)
    curyr = ((refdatet+cursecMF).strftime('%Y'))
    curmonth = ((refdatet+cursecMF).strftime('%m'))
    print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
    #initiating empty array for residual current calculation
    res_x = np.empty((len(xz),0))
    res_y = np.empty((len(xz),0))   

       
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        model_dimr.update()
        s1 = model_dfm.get_var('s1') # water level
        sa1 = model_dfm.get_var('sa1') # salinity
        # store the maximum water level per time step in column wise
        water_level = np.append(water_level, np.reshape(s1,(len(s1),1)), axis=1)
        salinity = np.append(salinity, np.reshape(sa1,(len(sa1),1)), axis=1)
        
        ## check to apply seedling establishment
        if curmonth == '01' and curyr_check == 0:
            print('prepare for seedlings establishment')
            res_x, res_y = collect_res( model_dfm, res_x, res_y)
            res_of_x = calculate_residual(res_x, model_dfm.get_time_step(), xz).reshape((len(res_x),1))
            res_of_y = calculate_residual(res_y, model_dfm.get_time_step(), xz).reshape((len(res_y),1))
            # create matrix (array)
            residual_is = np.hstack((res_of_x,res_of_y))     
        elif curyr_check != 0:
            if curmonth == '01' and curyr != curyr_check:
                print('prepare for seedlings establishment')
                res_x, res_y = collect_res( model_dfm, res_x, res_y)
                res_of_x = calculate_residual(res_x, model_dfm.get_time_step(), xz).reshape((len(res_x),1))
                res_of_y = calculate_residual(res_y, model_dfm.get_time_step(), xz).reshape((len(res_y),1))
                # create matrix (array)
                residual_is = np.hstack((res_of_x,res_of_y))       
            else:
                print(curyr, '/', curmonth, 'no seedlings establishment')
            
        # calculate the drag coefficient
        drag_coeff = newCalcDraginLoop(xyzw_cell_number, xyzw_nodes, xk, yk, 
                                       read_data, model_dfm, index_veg_cel)
        # drag_coeff = np.append(drag_coeff, addition)*ind 
        drag_coeff = np.append(drag_coeff, addition)
        # update the variable with new value
        model_dfm.set_var('Cdvegsp',drag_coeff)
        dts = model_dfm.get_time_step()
        t=t+dts
        print('Coupling ',str(ntime+1), 'run ', t, '/', coupling_period_model)
    
   
    med_sal = np.median(salinity, axis=1)
    print('Calculate median salinity in coupling',str(ntime+1))
    
    # 1.3. Prepare pandas for list of the seedling_finalpos
    
    # check condition
    if curmonth == '01' and curyr_check == 0:
        print('calculate seedlings position for year', curyr )
        seedling_finalpos = seedling_dispersal(xyzw_cell_number, index_veg_cel, 
                                               xyzw_nodes, xk, yk, read_data, 
                                               med_sal, residual_is, model_dfm, 
                                               reftime, cursec)
        # check seedlngs duplicate
        bool_series = seedling_finalpos.duplicated(keep='first')
        seedling_finalpos = seedling_finalpos[~bool_series]
        
        seedling_finalpos.to_csv(os.path.join(seedlings_out, 
                                'Seedling_Coupling_'+str(ntime+1)+'.txt'), 
                                 sep=',', index=False, header=True)
        # This is occupied for next looping checking
        current_sec = datetime.timedelta(seconds=model_dfm.get_current_time())
        # current simulation time multiply by MF
        current_secMF = current_sec*MorFac
        # current month based on simulation time (with MF)
        curyr_check = ((refdatet+current_secMF).strftime('%Y'))

    elif curyr_check != 0:
        if curmonth == '01' and curyr != curyr_check:
            print('calculate seedlings position for year', curyr )
            seedling_finalpos = seedling_dispersal(xyzw_cell_number, index_veg_cel, 
                                                   xyzw_nodes, xk, yk, read_data, 
                                                   med_sal, residual_is, model_dfm, 
                                                   reftime, cursec)
            # check seedlngs duplicate: https://www.machinelearningplus.com/pandas/pandas-duplicated/
            bool_series = seedling_finalpos.duplicated(keep='first')
            seedling_finalpos = seedling_finalpos[~bool_series]
            
            seedling_finalpos.to_csv(os.path.join(seedlings_out, 
                                    'Seedling_Coupling_'+str(ntime+1)+'.txt'), 
                                     sep=',', index=False, header=True)
            # This is occupied for next looping checking
            current_sec = datetime.timedelta(seconds=model_dfm.get_current_time())
            # current simulation time multiply by MF
            current_secMF = current_sec*MorFac
            # current month based on simulation time (with MF)
            curyr_check = ((refdatet+current_secMF).strftime('%Y'))
        else:
            print(curyr, '/', curmonth, 'no seedlings establishment')
    
    
    
    ### 2. convert the water_level from each time_step to each day
    # calculate the x axis of the array of the default run
    wl_shape = water_level.shape
    per_column = model_dfm.get_time_step()*MorFac
    time_linspace = np.linspace(per_column, coupling_period, wl_shape[1])
    # get the time step per column
    value_floor = np.floor(per_column/3600) # get the smaller step to get denser array (in hour)
    # check whether the value is more than 6 that makes it difficult to have an even array
    if (value_floor % 2) == 0:
        value_floor = value_floor
    elif(value_floor % 2) == 1 and value_floor > 6:
        value_floor = 6
    else:
        value_floor = value_floor

    # create an even x-axis    
    value_interp = int(coupling_period/3600/value_floor)
    time_interp = np.linspace(per_column, coupling_period, value_interp)

    # use interp1d to calculate the interpolated water level
    water_level_interp = np.empty((0, value_interp))
    for row in range(int(wl_shape[0])):
        f = interp1d(time_linspace,water_level[row,:])
        wl = f(time_interp)
        water_level_interp = np.append(water_level_interp, np.reshape(wl,(1,value_interp)), axis=0)

    col_num = int(24/value_floor) # equal to how many columns represent 1 day
    col_lookup = int(value_interp/col_num)
    # find daily maximum water level
    # wl_shape = water_level.shape
    h_wl = np.empty((0, col_lookup))
    # real time in hour is coupling_period/3600
    bb = np.empty((0, col_lookup))

    for ii in range(int(wl_shape[0])):
        cc = water_level_interp[ii,:]
        bb = []
        for aa in range(col_lookup):
            # bb_col = np.amax(cc[:,aa*col_num:aa*col_num+col_num])
            bb_col = np.amax(cc[aa*col_num:aa*col_num+col_num])
            bb.append(bb_col)
            # np.concatenate((bb,bb_col))
        bb = np.array(bb)
        bb = bb.reshape(1,col_lookup)
        h_wl = np.append(h_wl,bb,axis=0)
    
    ### 3. Calculate the probability based on the WoO and store the value in each cell
    # find median value of the h_wl for each cell number
    med_h_wl = np.median(h_wl, axis=1) 
    # find median value of salinity for each cell number
    med_sal = med_sal
    
    # calculate WoO probability value from h_wl (daily max water level)
    surv_val = np.empty(len(med_h_wl)) #initiate an empty array
    for ii in range(h_wl.shape[0]):
        fromheightcalc, Pvaluecalc = calcWOO(h_wl[ii,:],woo_inun) # get correlation of elevation and Probability
        surv_val[ii] = np.interp(med_h_wl[ii],fromheightcalc,Pvaluecalc)
        
    ## Filter seedlings based on the surv_val
    seedling_finalpos_2 = seedling_prob(seedling_finalpos, xyzw_cell_number, 
                                     xyzw_nodes, xk, yk, surv_val)
    seedling_finalpos = seedling_finalpos_2
    
    ### 4. Convert from data point to raster environment 
    # 4.1 Create raster from the surv-val
    surv_val_raster = np.column_stack((xz,yz,surv_val))
    concave_path = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
    if not os.path.exists(concave_path):
        os.makedirs(concave_path)
    surv_val_name = 'CH_surv_val'
    csv2ClippedRaster(concave_path, surv_val_raster, surv_val_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
    #return to home
    os.chdir(PROJ_HOME)
    
    #4.2 Create raster from med_sal
    sal_raster = np.column_stack((xz,yz,med_sal))
    sal_val_name = 'CH_sal_val'
    csv2ClippedRaster(concave_path, sal_raster, sal_val_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
    #return to home
    os.chdir(PROJ_HOME)
    
    
    ### 5. Tile the raster based on the tile in initialization folder
    # 5.1 Tile surv_val raster
    filename_surv = "tile_surv_"
    ras_clip = os.path.join(concave_path, surv_val_name+affix+'.tif')
    Sald3dNewRaster2Tiles(ras_clip, concave_path, tile_size_x, tile_size_y, dir_out, filename_surv)

    #5.2 Tile sal_val raster
    filename_sal = "tile_sal_"
    ras_clip = os.path.join(concave_path, sal_val_name+affix+'.tif')
    Sald3dNewRaster2Tiles(ras_clip, concave_path, tile_size_x, tile_size_y, dir_out, filename_sal)    
    
    ### 6. Create Mangrove Trees shp and tile
    # 6.1. createPointSHP(read_data, age_coupling0, concave_path, EPSG_Project) 
    createPointSHP(read_data, age_coupling, concave_path, EPSG_Project) 
    
    # 6.2. Tile the Master Trees and Save as XLS Input File
    shp_source = os.path.join(concave_path, 'coupling'+str(ntime+1)+'.shp') # location of the source tree and the master tree shp
    folder_loc = dir_out # location of the master tiles
    file_tile = os.path.join(dir_out,'tile_*.shp' )
    save_tiled_trees = os.path.join(MFON_Trees,'coupling'+str(ntime+1)) # location to save the tiled shp
    if not os.path.exists(save_tiled_trees):
        os.makedirs(save_tiled_trees)
      
    # clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name, a0, b0, a137, b137)
    New_clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name)
        
    ### 8. Create the env raster and tile raster files
    save_tiled_env = os.path.join(MFON_Env, 'coupling'+str(ntime+1))
    if not os.path.exists(save_tiled_env):
        os.makedirs(save_tiled_env)
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    # calc_surv = '"0*(A==-999)"'
    # calc_sal = '"0*(A>-999)+35"'
    # val_no_data_surv = 0
    # val_no_data_sal = 60
    
    # 8.1 Prepare Surv_Val Env Raster for MesoFON 
    filename_surv_env = filename_surv+'*'     
    end_name_surv = '_surv_'      
    SalNew_func_createRaster4MesoFON(concave_path,save_tiled_env, 
                                     filename_surv_env, dir_out, end_name_surv)
    
    # 8.2 Prepare Sal_Val Env Raster for MesoFON 
    filename_sal_env = filename_sal+'*' 
    end_name_sal = '_sal_'
    val_no_data_sal = 60
    SalNew_Sal_func_createRaster4MesoFON(concave_path,save_tiled_env, 
                                     filename_sal_env, dir_out, end_name_sal,
                                     val_no_data_sal)
    
    ### 9. Replace the content in Unrolled_Param, batch_params.xml, and parameters.xml
    # function to replace the content of the file
    # for filepatt in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
        
    Surv_Source, Sal_Source, Excel_Source = modifyParamMesoFON(MFON_HOME, MFON_Exchange, 
                                                save_tiled_env, save_tiled_trees, 
                                                Surv_Source, Sal_Source, 
                                                Excel_Source, ntime)
    
    ### 10. Run the MesoFON
    for filepatt in glob.iglob(os.path.join(MFON_HOME, 'tile_*')):
        # delete the existing instance1 in in the folder to prevent symlink errorp prior running
        try:
            send2trash.send2trash(os.path.relpath(os.path.join(filepatt,'instance_1'),PROJ_HOME))
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
        try:
            df = pd.read_csv(nama_a)
        except pd.errors.EmptyDataError:
            print( nama_a, " is empty")
        all_df.append(df)

    Concat_table = pd.concat(all_df)
    
    # 12.1. drop tick 0 year, only take 0.25
    Concat_table = Concat_table[Concat_table.tick > 0]
    run_is = 'Coupling_'+str(ntime+1) # change this with the real name
    
    # 12.2. use spatial in scipy to match the x,y of the mangroves and the age information.
    master_trees = gpd.read_file(os.path.join(concave_path+str('\\')+Path(concave_path).stem+'.shp'))
    age_coupling = calcAgeCoupling0(Concat_table, master_trees) # prepare the age_coupling for the Master Trees of Each Coupling procedure
    Concat_table['Age'] = age_coupling #update age after MesoFON run
    Concat_table.drop(['tick'], inplace=True, axis=1)
    Concat_table = Concat_table.reset_index(drop=True) # reset the index
    
    # 12.3. Concatenated table is saved as txt file
    Concat_table.to_csv(os.path.join(MFON_OUT_compile, run_is+'.txt'), sep=',', index=False, header=True)
    
    ### 13. Read the compile txt file and create the Cd
    read_data = Concat_table  
    
    print('End of coupling',str(ntime+1))
    print('Date now with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 

    
### End Loop
#Finalize the running
model_dimr.finalize()   




















# %%
