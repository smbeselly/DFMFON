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

#%% Reset all variable
# Don't forget to %reset %clear

#%% Input Folders and Files
PROJ_HOME = r'C:\Users\sbeselly\Downloads\02_1_fullmsl_sn_ext_200x200'
SYS_APP = r'C:\Users\sbeselly\Downloads\02_1_fullmsl_sn_ext_200x200\FnD3D'
D3D_HOME = r'C:\Users\sbeselly\Downloads\MASTER\2sebrian_20220518' #sal_veg-OK
gdal_loc = r'C:\Users\sbeselly\AppData\Local\miniconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbeselly\Downloads\MASTER\JDK11\bin\java.exe'
Mangr_SHP = '3m_04_full_mhwn-msl_200x200.shp'

D3D_Model = 'd3d_wave_r8_synth'
D3D_Domain = '20x410_net.nc' 
config_xml = 'r8_withwave_sn_ext.xml'
mdu_file = 'FlowFM_with_wave.mdu'                                             

#%% Import the packages and set the file
import numpy as np
import numpy.ma as ma
import bmi.wrapper                             
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
from dfm_tools.io.mdu import read_deltares_ini
from d3d_meso_mangro import createPointSHPmulti, n_calcAgeCoupling0   
from d3d_meso_mangro import modifyParamMesoFON 
from d3d_meso_mangro import csv2ClippedRaster, Sald3dNewRaster2Tiles
from d3d_meso_mangro import SalNew_func_createRaster4MesoFON 
from d3d_meso_mangro import New_clipSHPcreateXLSfromGPD, SalNew_Sal_func_createRaster4MesoFON
from d3d_meso_mangro import create_xzyzCellNumber, n_definePropVeg 
from d3d_meso_mangro import calcLevelCellCdveg, list_subsetCdveg
from d3d_meso_mangro import fill_elev_base, fill_elev_sdlg, add_cell_number
from d3d_meso_mangro import n_n_prepCalcDraginVect, n_n_calcDraginVect_fromPrep, calc_woo_avc_rzp_hs

from d3d_mangro_seeds import index_veg_cdveg
from d3d_mangro_seeds import seedling_dispersal, calculate_residual, collect_res
from d3d_mangro_seeds import elim_seeds_ced, seedling_prob_vect
 

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
from dateutil.rrule import rrule, DAILY #, MONTHLY

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
seedlings_before_woo = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings_before_woo')
if not os.path.exists(seedlings_before_woo):
    os.makedirs(seedlings_before_woo)
ddr = os.path.join(MFON_OUT,'Drag_Compile')
if not os.path.exists(ddr):
    os.makedirs(ddr)
#%% Settings
EPSG_Project = 3395 #EPSG worldwide mercator
couplingperiod = 90 #days
MorFac = 90 
woo_inun_av = 3 # inundation free period (days)
woo_inun_rh = 5
LLWL = -1.1
species_name = ['Avicennia_marina', 'Rhizopora_apiculata'] # must as list

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
ugrid_all = get_netdata(file_nc=grid_file)

xzyz_cell_number = create_xzyzCellNumber(xz, yz, model_dfm, ugrid_all)

#### Read Master Trees
master_trees = gpd.read_file(os.path.join(MFON_Trees, 'Master-Trees', Mangr_SHP))


#%% Read the compiled tree from the MesoFON Initialization Run and calculate drag_coefficient

MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')
## for restart
read_data = pd.read_csv(os.path.join(MFON_OUT_compile,'Coupling_0.txt')) #To restart from coupling 57
read_data['elev_base'] = np.nan #to prepare elev base column

## normal run
# use spatial in scipy to match the x,y of the mangroves and the age information.

age_coupling = read_data['Age']

# For loop for all of the cell number
index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
                                                             
# assume that the boundary flow nodes are located at the end of array
addition = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + 0.005  
addition_veg = copy.copy(addition)*0

drag_coeff = np.ones(model_dfm.get_var('ndx'))*0.005

# initiate calculate subsurface contribution of mangrove roots
bl_val = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
addition_bl = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_val[-1]
# model_dfm.get_var('bl')[bl_val.shape[0]-1:-1]
bl_val = np.append(bl_val, addition_bl)
# get base elevation info
fill_elev_base(model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)

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
coupling_period = couplingperiod*24*3600 
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
## force coupling_ntimeUse to be equal as in the DFMFON run
coupling_ntimeUse = 81 # force this to 20 years, change later
        
for ntime in range(int(coupling_ntimeUse)):
    # do the calculation for each coupling_ntime
    timetest = datetime.datetime.now() 
    print('Start the coupling',str(ntime+1),'computation')
    ### 1.1. run the DFM all the simulation time within ntime
    # find cells that have vegetation
    index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
    # predefine the pandas dataframe of the mangrove positions
    list_read_subset = list_subsetCdveg(ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
 
    # update the variable with new value taken from previous drag calculation
    model_dfm.set_var('Cdvegsp',drag_coeff)
    water_depth = np.empty((len(xz),0)) 
    salinity = np.empty((len(xz),0)) 
    bed_level = np.empty((len(xz),0))  
    ddrag = np.empty((len(xz),0))  ## to store the drag coeficient
    #https://www.delftstack.com/howto/numpy/python-numpy-empty-array-append/
    
    # calculate mangroves' subsurface contribution
    try:
        bl_model = model_dfm.get_var('bl') + bl_delta # the correct script
        model_dfm.set_var('bl',bl_model)
    except:
        print('no bed level update, model initiation')
        
    try:
        list_seed2sapl.append(seedling_finalpos)
        seedling_finalpos = pd.DataFrame() # to make sure no seedlings are fed when it is not the season
        seed2sapl = pd.concat(list_seed2sapl, axis=0)
        # add age to the appended pandas
        seed2sapl['Age'] = seed2sapl['Age']+add_seeds_age
        seed2sapl.to_csv(os.path.join(seedlngs_w_age_out, 
                                'Seedling_with_age_Coupling_'+str(ntime+1)+'.txt'), 
                                  sep=',', index=False, header=True)
        
        list_seed2sapl = []
        list_seed2sapl.append(seed2sapl)
        
        check_as_saplings = seed2sapl[seed2sapl['Age'] >= datetime.timedelta(days = 730)]

        if check_as_saplings.size > 0:
            print('Transfer seedlings as saplings at coupling', str(ntime+1))
            # select seedlings that have been transformed to saplings
            now_as_saplings = check_as_saplings.copy()
            now_as_saplings.loc[:, ['Age']] = (now_as_saplings['Age']/datetime.timedelta(days=365))
            now_as_saplings['dbh_cm'] = 1
            now_as_saplings['Height_cm'] = np.where(now_as_saplings['Species'] == 'Avicennia_marina',
                                                    141.80228, 144.72204)
            
            try:
                now_as_saplings.drop(['ParentPosX', 'ParentPosY','row',
                                      'Created Time Stamp', 'CrownSurfaceArea_m2',
                                      'cell_id', 'surv_val_av', 'surv_val_rh', 'cell_number'], 
                                     inplace=True, axis=1)
                print('saplings are Avicennia spp and Rhizopora spp')
            except:
                print('either Avicennia spp or Rhizopora spp only')            
            
                try:
                    now_as_saplings.drop(['ParentPosX', 'ParentPosY','row',
                                        'Created Time Stamp', 'CrownSurfaceArea_m2',
                                        'cell_id', 'surv_val_av', 'cell_number'], 
                                        inplace=True, axis=1)
                    print('domain contains Avicennia spp only')
                except:
                    now_as_saplings.drop(['ParentPosX', 'ParentPosY','row',
                                        'Created Time Stamp', 'CrownSurfaceArea_m2',
                                        'cell_id', 'surv_val_rh', 'cell_number'], 
                                        inplace=True, axis=1)
                    print('domain contains Rhizopora spp only')
            
            now_as_saplings['elev_base'] = np.nan
            index_spl = index_veg_cdveg(xzyz_cell_number, ugrid_all, now_as_saplings)
            fill_elev_base(model_dfm, ugrid_all, xzyz_cell_number, index_spl, now_as_saplings)
            
            # update the read data with new saplings as mature mangroves
            read_data = pd.concat([read_data, now_as_saplings], axis=0, ignore_index=True)
                        
            # reset list to use the filterd list from current selection
            # filter for less than 730
            filter_seeds2apl = seed2sapl[seed2sapl['Age'] <= datetime.timedelta(days = 730)]
            # list_seed2sapl = []
            seed2sapl = filter_seeds2apl.copy()
            list_seed2sapl = []
            list_seed2sapl.append(seed2sapl)
            
            len_rz = len(now_as_saplings[now_as_saplings['Species']=='Rhizopora_apiculata'])
            len_av = len(now_as_saplings[now_as_saplings['Species']=='Avicennia_marina'])
            print('Total number of saplings Avicennia: {} Rhizopora: {}:'.format(len_av, len_rz))
            
            now_as_saplings.to_csv(os.path.join(saplings_out, 
                                    'Saplings_Coupling_{}'.format(ntime+1)+'.txt'), 
                                      sep=',', index=False, header=True)
            
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
    
    ## test new January checking
    cur_date = refdatet+cursecMF
    nxt_secMF = cursecMF + datetime.timedelta(seconds=(coupling_period_model*MorFac)) #get date after coupling
    nxt_date = refdatet + nxt_secMF
    
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
    with_cell_num = add_cell_number(ugrid_all, xzyz_cell_number, read_data)
    avicennia_pr, RhDFrm_pr, vol_func, area_func = n_n_prepCalcDraginVect(with_cell_num, SYS_APP)
    rnveg_coeff, diaveg_coeff, stemheight_coeff = n_definePropVeg(avicennia_pr, RhDFrm_pr, vol_func, model_dfm)
    # copy new name reserved for tree detecting to adjust age on the MFON Trees
    mature_sapl = with_cell_num.copy()
    
    model_dfm.set_var('rnveg',rnveg_coeff)
    model_dfm.set_var('diaveg',diaveg_coeff)
    model_dfm.set_var('stemheight',stemheight_coeff)
    
    timetestend = datetime.datetime.now() 
    timeneeded = timetestend - timetest
    print('runtime initial coupling {}'.format(timeneeded))
    
    timeisnow = datetime.datetime.now()   
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        model_dimr.update()
        hs = model_dfm.get_var('hs') # water depth
        sa1 = model_dfm.get_var('sa1') # salinity
                                                                    
        water_depth = np.append(water_depth, np.reshape(hs,(len(hs),1)), axis=1)
        salinity = np.append(salinity, np.reshape(sa1,(len(sa1),1)), axis=1)
        
        ## check to apply seedling establishment
        if curmonth == '01' and curyr_check == 0:
            print('prepare for seedlings establishment')
            res_x, res_y = collect_res( model_dfm, res_x, res_y)
            blv = model_dfm.get_var('bl') # get the bed level
            bed_level = np.append(bed_level, np.reshape(blv,(len(blv),1)), axis=1)
            set_prep_seeds = True                     
        
        elif curyr_check != 0:
            try:
                if day_counts[1] < 16 and curyr != curyr_check:
                    print(curyr, '/', curmonth, 'no seedlings establishment \ndoes not satisfy min num of days') ## indicates that January days are not enough
                elif day_counts[1] < 16 and curyr == curyr_check:
                    print(curyr, '/', curmonth, 'no seedlings establishment \ndoes not satisfy min num of days') ## indicates that seedlings have been established in previous coupling
                elif day_counts[1] >= 16:
                    if 1 in chk_seed_prod and curmonth == '01':
                        print('prepare for seedlings establishment')
                        res_x, res_y = collect_res( model_dfm, res_x, res_y)
                        blv = model_dfm.get_var('bl') # get the bed level
                        bed_level = np.append(bed_level, np.reshape(blv,(len(blv),1)), axis=1)
                        set_prep_seeds = True                     
            except:
                print(curyr, '/', curmonth, 'no seedlings establishment')
        
        dts = model_dfm.get_time_step()
        t=t+dts
        timeis_loop = datetime.datetime.now()
        if t % model_dfm.get_time_step() == 0:
            # calculate the drag coefficient
            drag_coeff = n_n_calcDraginVect_fromPrep(avicennia_pr, RhDFrm_pr, vol_func, area_func, model_dfm, drag_coeff)
                                                                               
                                                        
            # update the variable with new value
            model_dfm.set_var('Cdvegsp',drag_coeff)
        timeend_loop = datetime.datetime.now()
        print('Coupling ',str(ntime+1), 'run ', t, '/', coupling_period_model)
        print('Runtime is ', timeend_loop-timeis_loop)
        
        cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) #https://stackoverflow.com/questions/775049/how-do-i-convert-seconds-to-hours-minutes-and-seconds
        # current simulation time multiply by MF
        cursecMF = cursec*MorFac
        # current month based on simulation time (with MF)
        curyr = ((refdatet+cursecMF).strftime('%Y'))
        curmonth = ((refdatet+cursecMF).strftime('%m'))
        print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
        
        ddrag = np.append(ddrag, np.reshape(drag_coeff,(len(sa1),1)), axis=1)
        
    timeisend = datetime.datetime.now()
    print('Runtime', 'Coupling ',str(ntime+1), 'is', timeisend-timeisnow)
    curyr_check = curyr

    ddrag_med = np.median(ddrag, axis=1)
    np.savetxt(os.path.join(ddr, 'Coupling_{}.txt'.format(ntime+1)), ddrag_med, delimiter=',')
    
    if set_prep_seeds == True :
        res_of_x = calculate_residual(res_x, model_dfm.get_time_step(), xz).reshape((len(res_x),1))
        res_of_y = calculate_residual(res_y, model_dfm.get_time_step(), xz).reshape((len(res_y),1))
        # create matrix (array)
        residual_is = np.hstack((res_of_x,res_of_y))
        print('calculate residual current in coupling',str(ntime+1))
        # reset set_prep_seeds                     
        set_prep_seeds = False                      
    else:
        residual_is = []
    med_sal = np.median(salinity, axis=1)
    print('Calculate median salinity in coupling',str(ntime+1))
    
    # 1.3. Prepare pandas for list of the seedling_finalpos
    
    # check condition
    if len(residual_is) != 0:
        print('calculate seedlings position for year', curyr )
        seedling_finalpos = seedling_dispersal(xzyz_cell_number, index_veg_cel, ugrid_all, read_data, 
                               med_sal, residual_is, model_dfm, reftime, cursecMF)
        # check seedlngs duplicate
        bool_series = seedling_finalpos.duplicated(keep='first')
        seedling_finalpos = seedling_finalpos[~bool_series]
        # save seedling_finalpos before filter due to survival rate (WoO)
        seedling_finalpos.to_csv(os.path.join(seedlings_before_woo, 
                                'Seedling_before_Coupling_{}_{}'.format(ntime+1,curyr)+'.txt'), 
                                 sep=',', index=False, header=True)
    else:
        print(curyr, '/', curmonth, 'no seedlings establishment')
    
    
    ### 2. convert the water_level from each time_step to each day
    # calculate the x axis of the array of the default run
    def calc_h_wl(water_level, model_dfm, MorFac, coupling_period):
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
    
        return h_wl
    
    hs_wl = calc_h_wl(water_depth, model_dfm, MorFac, coupling_period)

    ### 3. Calculate the probability based on the WoO and store the value in each cell
    survval_hs = calc_woo_avc_rzp_hs(hs_wl, model_dfm)                                                                    
    surv_val_av = survval_hs[:,0]                                                                                                           
    surv_val_rh = survval_hs[:,1]
    
    # find median value of salinity for each cell number
    med_sal = med_sal
       
    if seedling_finalpos.size > 0:
        seedling_finalpos_2 = seedling_prob_vect(seedling_finalpos, xzyz_cell_number, 
                                            ugrid_all, surv_val_av, surv_val_rh)

        seedling_finalpos = seedling_finalpos_2
        try: # this is the control when all of the seedlings are died 
            # add base elevation for seedlings
            seedling_finalpos['elev_base'] = np.nan
            index_sdl = index_veg_cdveg(xzyz_cell_number, ugrid_all, seedling_finalpos)
            fill_elev_sdlg(bed_level[:,0], ugrid_all, xzyz_cell_number, index_sdl, seedling_finalpos)
            # filter seedlings based on critical erosion depth
            seedling_finalpos = elim_seeds_ced(bed_level, ugrid_all, xzyz_cell_number, index_sdl, seedling_finalpos)
        except:
            print('all of seedlings are died')
        seedling_finalpos.to_csv(os.path.join(seedlings_out, 
                                'Seedling_Coupling_{}_{}'.format(ntime+1,curyr)+'.txt'), 
                                 sep=',', index=False, header=True)

    
    ### 4. Convert from data point to raster environment 
    # 4.1 Create raster from the surv-val    
    surv_val_raster = np.column_stack((xz,yz,surv_val_av))
    concave_path = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
    if not os.path.exists(concave_path):
        os.makedirs(concave_path)
    surv_val_name = 'CH_surv_val'
    csv2ClippedRaster(concave_path, surv_val_raster, surv_val_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
    #return to home
    os.chdir(PROJ_HOME)

    ## raster for surv_val_rh
    surv_val_raster_rh = np.column_stack((xz,yz,surv_val_rh))
    concave_path_rh = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
    if not os.path.exists(concave_path_rh):
        os.makedirs(concave_path_rh)
    surv_val_name_rh = 'CH_surv_val_rh'
    csv2ClippedRaster(concave_path_rh, surv_val_raster_rh, surv_val_name_rh, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
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
    createPointSHPmulti(read_data, age_coupling, concave_path, EPSG_Project) 
    
    # 6.2. Tile the Master Trees and Save as XLS Input File
    shp_source = os.path.join(concave_path, 'coupling'+str(ntime+1)+'.shp') # location of the source tree and the master tree shp
    folder_loc = dir_out # location of the master tiles
    file_tile = os.path.join(dir_out,'tile_*.shp' )
    save_tiled_trees = os.path.join(MFON_Trees,'coupling'+str(ntime+1)) # location to save the tiled shp
    if not os.path.exists(save_tiled_trees):
        os.makedirs(save_tiled_trees)
      
                                                                                                        
    New_clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name)
        
    ### 8. Create the env raster and tile raster files
    save_tiled_env = os.path.join(MFON_Env, 'coupling'+str(ntime+1))
    if not os.path.exists(save_tiled_env):
        os.makedirs(save_tiled_env)
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    
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
        # only calculate MesoFON if trees exist 
        if Path(os.path.join(save_tiled_trees,Path(filepatt).stem[:-12]+'.shp')).is_file():  
            print(Path(filepatt).stem[:-6])    
            
            # delete the existing instance1 in in the folder to prevent symlink errorp prior running
            try:
                # send2trash.send2trash(os.path.relpath(os.path.join(filepatt,'instance_1'),PROJ_HOME))
                directory = os.path.join(filepatt,'instance_1')
                files_in_directory = os.listdir(directory)
                filtered_files = [file for file in files_in_directory if file.startswith("MF")]
                for file in filtered_files:
                    path_to_file = os.path.join(directory, file)
                    send2trash.send2trash(path_to_file)
            except OSError as e:
                print("Instance 1 is already deleted before this command: %s : %s" % (os.path.join(filepatt,'instance_1'), e.strerror))
													
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
        if nama != []:
            shutil.copyfile(nama[0], os.path.join(MFON_OUT_tile,Path(nama[0]).stem +' Coupling_'+str(ntime+1)+'.txt'))
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
            df = pd.DataFrame()
        all_df.append(df)

    Concat_table = pd.concat(all_df)
    
    # 12.1. drop tick 0 year, only take 0.25
    Concat_table = Concat_table[Concat_table.tick > 0]
    run_is = 'Coupling_'+str(ntime+1) # change this with the real name
    
    # 12.2. use spatial in scipy to match the x,y of the mangroves and the age information.
    # master_trees = gpd.read_file(os.path.join(concave_path+str('\\')+Path(concave_path).stem+'.shp'))
    # age_coupling = calcAgeCoupling0(Concat_table, master_trees) # prepare the age_coupling for the Master Trees of Each Coupling procedure
    age_coupling = n_calcAgeCoupling0(Concat_table, mature_sapl) 
    Concat_table['Age'] = age_coupling.to_numpy() #update age after MesoFON run
    Concat_table.drop(['tick'], inplace=True, axis=1)
    Concat_table = Concat_table.reset_index(drop=True) # reset the index
                                                                                           
    
    # 12.3. Concatenated table is saved as txt file
    Concat_table.to_csv(os.path.join(MFON_OUT_compile, run_is+'.txt'), sep=',', index=False, header=True)
    # save bottom depth as text to facilitate the ongoing simulation visualisation 
    np.savetxt(os.path.join(botdepth_out, run_is+'.txt'), model_dfm.get_var('bl'), delimiter=",")
    
    
    ### 13. Read the compile txt file to prepare for the next iteration
    read_data = Concat_table  
    
    # calculate below ground biomass contribution to bed level
    bl_mangro = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
    addition_bl_mangro = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_mangro[-1]
    bl_mangro = np.append(bl_mangro, addition_bl_mangro)
    # calculate the delta after every coupling 
    bl_delta = bl_mangro - bl_val
    # save current calcLevelCell value for the next iteration
    bl_val = bl_mangro.copy()
    
    timeislast = datetime.datetime.now()
    print('End of coupling',str(ntime+1))
    print('Date now with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
    print('Runtime', 'End of coupling ',str(ntime+1), 'is', timeislast-timeisnow)
    
    # 12.4 Read Canopy Out, Compile, and Save as txt
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
            shutil.copyfile(nama[0], 
                            os.path.join(MFON_OUT_tile,Path(nama[0]).stem +run_is+'.txt'))
            namac.append(nama[0])

    ### Compile the results to compile canopy folder
    MFON_OUT_compile_can = os.path.join(MFON_OUT,'Compile_Canopy')
    if not os.path.exists(MFON_OUT_compile_can):
        os.makedirs(MFON_OUT_compile_can)

    all_dfc = []    
    for nama_a in namac:
        df = pd.read_csv(nama_a)
        all_dfc.append(df)

    Concat_tablec = pd.concat(all_dfc)
    Concat_tablec = Concat_tablec.reset_index(drop=True) # reset the index
    # Concatenated table is saved as txt file
    Concat_tablec.to_csv(os.path.join(MFON_OUT_compile_can,run_is+'.txt'), 
                        sep=',', index=False, header=True)

   
### End Loop
#Finalize the running
model_dimr.finalize() 

