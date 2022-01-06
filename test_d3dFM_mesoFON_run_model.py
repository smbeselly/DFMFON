# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 20:00:56 2022

@author: sbe002

TODO belum nambah script untuk tambahan sediment karena biomass (daun jatuh, dll)
"""

#%% Import the necessary packages, set the file path, and input files

import numpy as np
import numpy.ma as ma
import bmi.wrapper
import ctypes
import cmocean.cm
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os
import datetime
import pandas as pd
from scipy import integrate
import faulthandler
faulthandler.enable()
import sys
# print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes, calcDragCoeff    

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

## Set the paths for dll-files and input-files for DFM
D3D_HOME = os.path.join(r'C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64')
MFON_HOME = os.path.join(r'Model-Execute/MesoFON')
D3D_workdir = os.path.join(r'Model-Execute/D3DFM/FunnelMorphMF30') # model funnel with morpho
MFON_Exchange = os.path.join(r'Model-Exchange')
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
if not os.path.exists(MFON_Env):
    os.makedirs(MFON_Env)

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

#%% Define the further path

## No need to update this,, except the config_file

dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm','bin')
dflowfm_engine = os.path.join(dflowfm_path, 'dflowfm.dll')
config_file = os.path.join(D3D_workdir, 'FunnelMorphMF30.xml') # funnel morpho

mdu_file = os.path.join(D3D_workdir, 'dflowfm', 'FlowFM.mdu')
grid_file = os.path.join(D3D_workdir, 'dflowfm', 'Delta_Schematized_funnel_net.nc')
figsavefolder= os.path.join(r'Model-Out/D3DFM/Figs')

# Uncomment this if you want to run with console
# print("D3D_HOME :" + D3D_HOME)
# print("dimr_path :" + dimr_path)
# print("config_file:" + config_file)
# print("MFON_HOME:" + MFON_HOME)
# print("mfon_path:" + mfon_path)
# print("mfon_config_file:" + mfon_config_file)

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

## Define DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
