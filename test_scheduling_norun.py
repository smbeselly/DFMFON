# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:06:42 2022

@author: sbe002
"""

#%% Input Folders and Files
PROJ_HOME = r'C:\Users\sbe002\Downloads\Research_Run\Model_A_test_scheduling'
# PROJ_HOME = r'D:\Git\d3d_meso'
SYS_APP = r'D:\Git\d3d_meso/FnD3D'
D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\2sebrian_20220518' #sal_veg-OK
gdal_loc = r'C:\Users\sbe002\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'

D3D_Model = 'd3d_meso_run6_scenario_A'
D3D_Domain = 'Grid_Funnel_1_net.nc'
config_xml = 'd3d_meso_run6.xml'
# D3D_Model = 'model_run_6_flat'
# D3D_Domain = 'Grid_Funnel_20_by_20_net.nc'
# config_xml = 'model_run_6_flat.xml'
mdu_file = 'FlowFM.mdu'

Mangr_SHP = 'Tip_Saplings_geser_2.shp'
# Mangr_SHP = 'Tip_Saplings_geser_2_arah_y.shp'

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
# Supress/hide the warning
np.seterr(invalid='ignore')

from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc import get_netdata
# from dfm_tools.io.polygon import Polygon
from dfm_tools.io.mdu import read_deltares_ini
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes #, calcDragCoeff 
from d3d_meso_mangro import calcWOO, calcAgeCoupling0, createPointSHP #, createXLSfromSHP  
from d3d_meso_mangro import modifyParamMesoFON #, createRaster4MesoFON, calcDragInLoop
from d3d_meso_mangro import csv2ClippedRaster, Sald3dNewRaster2Tiles #, clipSHPcreateXLSfromGPD
from d3d_meso_mangro import SalNew_func_createRaster4MesoFON, newCalcDraginLoop
from d3d_meso_mangro import New_clipSHPcreateXLSfromGPD, SalNew_Sal_func_createRaster4MesoFON
from d3d_meso_mangro import initCalcDraginLoop, list_subset, calcLevelCell
from d3d_meso_mangro import create_xzyzCellNumber, initCalcDraginLoopCdveg,definePropVeg
from d3d_meso_mangro import calcLevelCellCdveg, list_subsetCdveg,newCalcDraginLoopCdveg

from d3d_mangro_seeds import index_veg_cdveg, seedling_establishment
from d3d_mangro_seeds import seedling_dispersal, calculate_residual, collect_res
from d3d_mangro_seeds import seedling_prob, elim_seeds_surv

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
import copy
from dateutil.rrule import rrule, MONTHLY, DAILY

# to log the print to  : https://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("logfile_02_seedling.log", "a")
   
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
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM',D3D_Model) # model funnel with morpho
# MFON_JAR = os.path.join(MFON_HOME, 'complete_model.jar')
MFON_LocalBatchRunner = os.path.join(MFON_HOME,'local_batch_run.properties')
gdal_path = os.path.join(gdal_loc)
JAVAREP = os.path.join(JAVA_Exe)

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')
DFM_OUT = os.path.join(PROJ_HOME,'Model-Out','D3DFM')

dir_out = os.path.join(MFON_Exchange, 'Initialization')

figsavefolder= os.path.join(PROJ_HOME,'Model-Out','Figures')
if not os.path.exists(figsavefolder):
    os.makedirs(figsavefolder)
seedlings_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings')
if not os.path.exists(seedlings_out):
    os.makedirs(seedlings_out)
saplings_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Saplings')
if not os.path.exists(saplings_out):
    os.makedirs(saplings_out)
seedlngs_w_age_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings_with_age')
if not os.path.exists(seedlngs_w_age_out):
    os.makedirs(seedlngs_w_age_out)
botdepth_out = os.path.join(DFM_OUT, 'Bottom_Depth')
if not os.path.exists(botdepth_out):
    os.makedirs(botdepth_out)
#%% Settings
EPSG_Project = 32749 # EPSG code for WGS84/ UTM Zone 49S (Porong case study)
coupling_period = 90 #days, actually 90 days
MorFac = 30 
woo_inun = 3 # inundation free period (days)
LLWL = -1.2
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

#%% Get time reference from the model

## get the reference date and starttime of model
getmdu = read_deltares_ini(mdu_file)
refdate = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'RefDate')]['value']
tstart = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'TStart')]['value']
tend = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'TStop')]['value']

# parse dfm's time in string to datetime var in Python
refdatet = parser.parse(refdate.iloc[0])
tstartt = datetime.timedelta(seconds=float(tstart.iloc[0]))
tendd = datetime.timedelta(seconds=float(tend.iloc[0]))
#reference time is
reftime = refdatet+tstartt
timeend = refdatet+tendd
timendwmorf = refdatet+ (tendd*MorFac)

print('start simulation time', reftime)
print('end of simulation in DFM', timeend)
print('end of simulation with MorFac', timendwmorf)

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

#%% Fake DFM Runner
timing = np.arange(0, int(tend)+200,200)

# timesecs = np.empty(shape=[0,np.size(timing)])
timesecs = []

#%% run DFM test
for ntime in range(int(coupling_ntimeUse)):
# for ntime in range(int(6)):
    # do the calculation for each coupling_ntime
    print('Start the coupling',str(ntime+1),'computation')
    
    ### 1.2. Check for seedling establishment
    # seedling establishment only occur during fruiting season (January)
    # current simulation time is
    # cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) #https://stackoverflow.com/questions/775049/how-do-i-convert-seconds-to-hours-minutes-and-seconds
    try:
        curseco = datetime.timedelta(seconds=max(timesecs))
    except:
        curseco = datetime.timedelta(seconds=0)
    # current simulation time multiply by MF
    cursecMF = curseco*MorFac
    # current month based on simulation time (with MF)
    curyr = ((refdatet+cursecMF).strftime('%Y'))
    curmonth = ((refdatet+cursecMF).strftime('%m'))
    print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
    #initiating empty array for residual current calculation
    # res_x = np.empty((len(xz),0))
    # res_y = np.empty((len(xz),0))   
    
    ## test new January checking
    cur_date = refdatet+cursecMF
    nxt_secMF = cursecMF + datetime.timedelta(seconds=(coupling_period_model*MorFac)) #get date after coupling
    nxt_date = refdatet + nxt_secMF
    
    # chk_seed_prod = [dayis.month for dayis in rrule(MONTHLY, dtstart=cur_date, until=nxt_date)]
    chk_seed_prod = [dayis.month for dayis in rrule(DAILY, dtstart=cur_date, until=nxt_date)]
    chk_seed_prod = list(set(chk_seed_prod))
    ## Control for Duplicate Month of January (when it leaps)
    def chck_full_jan():
        the_list = list(rrule(DAILY, dtstart=cur_date, until=nxt_date))
        day_counts = {}
        for day in the_list:
            day_counts[day.month] = day_counts.get(day.month, 0)+1
        return day_counts    

    day_counts = chck_full_jan()    

    # timeisnow = datetime.datetime.now()   
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        # model_dimr.update()
        # s1 = model_dfm.get_var('s1') # water level
        # sa1 = model_dfm.get_var('sa1') # salinity
        # # store the maximum water level per time step in column wise
        # water_level = np.append(water_level, np.reshape(s1,(len(s1),1)), axis=1)
        # salinity = np.append(salinity, np.reshape(sa1,(len(sa1),1)), axis=1)
        
        ## check to apply seedling establishment
        if curmonth == '01' and curyr_check == 0:
            print('prepare for seedlings establishment')
            # res_x, res_y = collect_res( model_dfm, res_x, res_y)
            dataframe = pd.DataFrame(list())
            dataframe.to_csv(os.path.join(seedlings_out,
                             '{}_{}'.format(curyr, curmonth)+'_throw.txt'))
        
        elif curyr_check != 0:
            try:
                if day_counts[1] < 16 and curyr != curyr_check:
                    print(curyr, '/', curmonth, 'no seedlings establishment \ndoes not satisfy min num of days') ## indicates that January days are not enough
                elif day_counts[1] < 16 and curyr == curyr_check:
                    print(curyr, '/', curmonth, 'no seedlings establishment \ndoes not satisfy min num of days') ## indicates that seedlings have been established in previous coupling
                elif day_counts[1] >= 16:
                    if 1 in chk_seed_prod and curmonth == '01':
                        print('prepare for seedlings establishment')
                        # res_x, res_y = collect_res( model_dfm, res_x, res_y) 
                        dataframe = pd.DataFrame(list())
                        dataframe.to_csv(os.path.join(seedlings_out,
                                         '{}_{}'.format(curyr, curmonth)+'_throw.txt'))
                        
            except:
                print(curyr, '/', curmonth, 'no seedlings establishment')

     
        # elif curyr_check != 0:
        #     if curmonth == '01' and curyr != curyr_check:
        #         print('prepare for seedlings establishment')
        #         res_x, res_y = collect_res( model_dfm, res_x, res_y)
        #     else:
        #         print(curyr, '/', curmonth, 'no seedlings establishment')
        
        # dts = model_dfm.get_time_step()
        dts = 200
        t=t+dts
        
        # if t % model_dfm.get_time_step() == 0:
        #     # calculate the drag coefficient
        #     drag_coeff = newCalcDraginLoopCdveg(model_dfm,xzyz_cell_number, 
        #                                         index_veg_cel,list_read_subset)
        #     drag_coeff = np.append(drag_coeff, addition)
        #     # update the variable with new value
        #     model_dfm.set_var('Cdvegsp',drag_coeff)
        
        print('Coupling ',str(ntime+1), 'run ', t, '/', coupling_period_model)

        # cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) #https://stackoverflow.com/questions/775049/how-do-i-convert-seconds-to-hours-minutes-and-seconds
        # current simulation time multiply by MF
        cursec = datetime.timedelta(seconds=t+curseco.total_seconds())
        cursecMF = cursec*MorFac
        # current month based on simulation time (with MF)
        curyr = ((refdatet+cursecMF).strftime('%Y'))
        curmonth = ((refdatet+cursecMF).strftime('%m'))
        print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
    timesecs.append(curseco.total_seconds()+t)  
    # timeisend = datetime.datetime.now()
    # print('Runtime', 'Coupling ',str(ntime+1), 'is', timeisend-timeisnow)
    curyr_check = curyr
   

    