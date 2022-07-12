## the source of the vegetation is coupling_0
# an exagerated vegetation value (and abnormal)

#%% Input Folders and Files
PROJ_HOME = r'D:\Git\d3d_meso'
SYS_APP = r'D:\Git\d3d_meso/FnD3D'
D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\2sebrian_20220518' #sal_veg-OK
gdal_loc = r'C:\Users\sbe002\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'

D3D_Model = 'model_run_6_flat'
D3D_Domain = 'Grid_Funnel_20_by_20_net.nc'
config_xml = 'model_run_6_flat.xml'
mdu_file = 'FlowFM.mdu'

# Mangr_SHP = 'geserMangroveAgeMerged.shp'
Mangr_SHP = 'Tip_Saplings_geser_2_arah_y.shp'

## this is additional for exagerated veg
is_veg = r"D:/Git/d3d_meso/00_for_test_coup0.txt"

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
botdepth_out = os.path.join(DFM_OUT, 'Bottom_Depth')
if not os.path.exists(botdepth_out):
    os.makedirs(botdepth_out)

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
# xyzw_cell_number = create_xyzwCellNumber(xzw,yzw,mesh_face_x,mesh_face_y)
# xyzw_nodes = create_xyzwNodes(mesh_face_nodes,xyzw_cell_number)
xzyz_cell_number = create_xzyzCellNumber(xz, yz, model_dfm, ugrid_all)

#### Read Master Trees
master_trees = gpd.read_file(os.path.join(MFON_Trees, 'Master-Trees', Mangr_SHP))

#### Read the polygon pli and add the indices
# pli = Polygon.fromfile(os.path.join(D3D_workdir,'dflowfm','vege.pli'))
# path = mpltPath.Path(pli[0][0])
# ind = path.contains_points(np.transpose((xz,yz))).astype(int)

# Read the compiled tree from the MesoFON Initialization Run and calculate drag_coefficient

MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')

read_data = pd.read_csv(is_veg)

#% Get time reference from the model

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

#% Loop the Coupling
coupling_period_model = 43200

#%% assign the Cd veg test one value but position is from read_data
# find cells that have vegetation

### This script has realistic drag but static rnveg, diamveg, and stemheight

# test if we want to shift the veg position
read_data['GeoRefPosX'] = read_data['GeoRefPosX']+380
##################
ind_is = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
addition_ind = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi')))
ind_is = np.append(ind_is, addition_ind)
######################
index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
drag_coeff = initCalcDraginLoopCdveg(xzyz_cell_number, model_dfm, ugrid_all, 
                                     index_veg_cel, read_data)
# assume that the boundary flow nodes are located at the end of array
addition = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + 0.005  
addition_veg = copy.copy(addition)*0
# drag_coeff = np.append(drag_coeff, addition)*ind
drag_coeff = np.append(drag_coeff, addition)
################
## initialisation of vegetation variables
rnveg=model_dfm.get_var('rnveg')  # dichtheid
diaveg=model_dfm.get_var('diaveg') # diameter
stemheight=model_dfm.get_var('stemheight') # hoogte
ndx = model_dfm.get_var_shape('s1')

rnveg=np.full(ndx,1.0)
diaveg=np.full(ndx,0.2)
stemheight=np.full(ndx,7.0)

model_dfm.set_var('rnveg',ind_is*rnveg)
model_dfm.set_var('diaveg',ind_is*diaveg)
model_dfm.set_var('stemheight',ind_is*stemheight)
model_dfm.set_var('Cdvegsp',drag_coeff)
#%% run model ## write as output_veg_static

t=0

while t<coupling_period_model:
    model_dfm.update()
    dts = model_dfm.get_time_step()
    t=t+dts
    
#Finalize the running
model_dfm.finalize() 


#%% assign the Cd veg test from read_data value and compare without veg

### This script has realistic drag, rnveg, diamveg, and stemheight

# test if we want to shift the veg position
read_data['GeoRefPosX'] = read_data['GeoRefPosX']+380

######################
drag_coeff = initCalcDraginLoopCdveg(xzyz_cell_number, model_dfm, ugrid_all, 
                                     index_veg_cel, read_data)
# assume that the boundary flow nodes are located at the end of array
addition = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + 0.005  
addition_veg = copy.copy(addition)*0
# drag_coeff = np.append(drag_coeff, addition)*ind
drag_coeff = np.append(drag_coeff, addition)
################

# find cells that have vegetation
index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
# predefine the pandas dataframe of the mangrove positions
list_read_subset = list_subsetCdveg(ugrid_all, xzyz_cell_number, index_veg_cel, read_data)

rnveg_coeff, diaveg_coeff, stemheight_coeff = definePropVeg(xzyz_cell_number, 
                        model_dfm, ugrid_all, index_veg_cel, read_data, addition_veg)

################

model_dfm.set_var('rnveg',rnveg_coeff)
model_dfm.set_var('diaveg',diaveg_coeff)
model_dfm.set_var('stemheight',stemheight_coeff)

# update the variable with new value taken from previous drag calculation
model_dfm.set_var('Cdvegsp',drag_coeff)

#%% run model ## write as output_veg_dynamic

t=0

while t<coupling_period_model:
    model_dfm.update()
    dts = model_dfm.get_time_step()
    t=t+dts
    
#Finalize the running
model_dfm.finalize() 
