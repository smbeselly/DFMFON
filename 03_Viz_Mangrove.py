# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 13:19:01 2022

This is the script to access the simulation results

@author: sbe002
"""
#%% Input Folders and Files
PROJ_HOME = r'E:\d3d_meso_run5_scenario_A_pcsleep_at_run_24sep'
SYS_APP = r'E:\d3d_meso_run5_scenario_A_pcsleep_at_run_24sep/FnD3D'
D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\oss_artifacts_x64_140691'
gdal_loc = r'D:\Program_Files\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'

D3D_Model = 'restart_d3d_meso_run5_scenario_A'
D3D_Domain = 'Grid_Funnel_1_net.nc'
config_xml = 'd3d_meso_run5_scenario_A.xml'
mdu_file = 'FlowFM.mdu'

Mangr_SHP = 'Saplings_run5A_Density_0.1.shp'

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
# from dfm_tools.io.polygon import Polygon
from dfm_tools.io.mdu import read_deltares_ini
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes #, calcDragCoeff 
from d3d_meso_mangro import calcWOO, calcAgeCoupling0, createPointSHP #, createXLSfromSHP  
from d3d_meso_mangro import modifyParamMesoFON #, createRaster4MesoFON, calcDragInLoop
from d3d_meso_mangro import csv2ClippedRaster, Sald3dNewRaster2Tiles #, clipSHPcreateXLSfromGPD
from d3d_meso_mangro import SalNew_func_createRaster4MesoFON, newCalcDraginLoop
from d3d_meso_mangro import New_clipSHPcreateXLSfromGPD, SalNew_Sal_func_createRaster4MesoFON
from d3d_meso_mangro import initCalcDraginLoop, list_subset
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
import natsort
import seaborn as sns
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
seedlings_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings')

shp_clip = os.path.join(dir_out, 'CH_shp')
concave_name = 'CH_bathy_'
affix = '_clipped'

#%% Access the results

MFON_Compile = os.path.join(MFON_OUT, 'Compile')
MFON_Seedlings = os.path.join(MFON_OUT, 'Seedlings')

## acquiring data from Compile folder
# get the list of txt file
listcomp = []
for filecomp in glob.iglob(os.path.join(MFON_Compile,'*.txt')):
    # print(Path(filecomp).stem)
    # ad = pd.read_csv(filecomp)
    listcomp.append(Path(filecomp).stem)
    # listcomp.append(ad)

# sort the text file to ascending mode
listcompt_sorted = natsort.natsorted(listcomp)

# read the compile txt as pandas and put as list
listcomp_baca =[]
for list_is in listcompt_sorted:
    ad = pd.read_csv(os.path.join(MFON_Compile,list_is+'.txt'))
    listcomp_baca.append(ad)
    
#%% visualisation
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
import cmocean

file_nc_map = os.path.join(D3D_workdir, 'dflowfm', 'output', mdu_file[:-4]+'_map.nc')
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_map)
times_pd = get_timesfromnc(file_nc=file_nc_map)

#plot net/grid
ugrid_all = get_netdata(file_nc=file_nc_map)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_aspect('equal')

#plot water level on map
# data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_s1', timestep=3)#, multipart=False)
data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_s1', timestep=3)#, multipart=False)
clmap = cmocean.cm.delta
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_wl[0,:], ax=None, linewidth=0.5, cmap=clmap)
# pc.set_clim([-0.5,1])
fig.colorbar(pc, ax=ax)
# ax.set_title('%s (%s)'%(data_frommap_wl.var_varname, data_frommap_wl.var_ncvarobject.units))
ax.set_aspect('equal')

# test seaborn
sns.set_theme(style="dark")
pakai = listcomp_baca[3]
sns.relplot(x='GeoRefPosX', y='GeoRefPosY', hue='Species' ,size='Height_cm', 
            sizes=(40, 400), alpha=0.5, palette='muted', data=pakai)
