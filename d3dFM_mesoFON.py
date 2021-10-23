# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Script to Coupling Delft3D and MesoFON

#%%
# Import the necessary packages
import numpy as np
import bmi.wrapper
import ctypes
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os
import datetime
from scipy import integrate
import faulthandler
faulthandler.enable()

#%%
## Initiate the BMI Delft3D
## Paths for dll-files and input-files for DFM
D3D_HOME = os.path.join(r'D:\Delft3D FM\FM_Compiled_1.2.8.62485\x64')
MFON_HOME = os.path.join(r'D:\Program_Files\eclipse-workspace\meso_FON')
workdir = os.path.join(r'D:\A_PhD\03. Bahan PhD\07. Olah Data\17._BMI_DFM')

# engine and config setting DFM
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm', 'bin', 'dflowfm.dll')
config_file = os.path.join(workdir, 'dimr_config.xml')

mdu_file = os.path.join(workdir, 'fm', 'sample.mdu')
grid_file = os.path.join(workdir, 'fm', 'sample_net.nc')
figsavefolder= os.path.join(r'D:\A_PhD\03. Bahan PhD\07. Olah Data\17._BMI_DFM\figs')

# engine and config setting MesoFON
# TODO
# mfon_path = os.path.join(MFON_HOME, 'src', 'meso_FON', 'DFM_Wrapper.java') ## maybe not necessary anymore
# mfon_config_file = os.path.join(MFON_HOME, 'batch', 'batch_params.xml')

print("D3D_HOME :" + D3D_HOME)
print("dimr_path :" + dimr_path)
print("config_file:" + config_file)
print("MFON_HOME:" + MFON_HOME)
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
model_dfm = bmi.wrapper.BMIWrapper(engine=dflowfm_path, configfile=mdu_file)

## Define and initialise DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()
print ('DFM initialized')

# ## Define MesoFON Wrapper
# model_mfon = bmi.wrapper.BMIWrapper(engine=mfon_path, config_file=mfon_config_file)
# model_mfon.initialize()
# print('MesoFON initialized')
#%%
## Read the D3D Domain and Map it to create the 'world' for MesoFON

## Get the limit of the approximated LLWL and build the outer line
# This is as the base for 'salinity approach' to kill the trees that are 
# located outside the allowed area

# =============================================================================
# # ~ try splitting the shp files and the raster
# # ~ create function based on this recipe "splitD3Dmap"
# Already done in testSplitD3D
# It creates a tiles of raster file with specific interval in meters or pixel
# TODO
# Change this as a function
# def splitD3Dmap():
#     pass
# =============================================================================

#%%
### Loop

    ## Read the trees location and create the polygon for Trachytope
    
    ## Do the first run of Delft3D
    
    # Pause the simulation to a certain delta t coupling (3 months)
    
    ## Script to process the output with BMI? or other default Python processing
    ## package for DFM
    
    ## Retrieve the D3DFM map based on WoO
    # run the "splitD3Dmap"
        # It will clip the map based on the master shapefile
        # Create raster based on the WoO as a Salinity, on the masterClipped 
        # However, it will only mark the unhabitable area as extreme salinity
        # ----
        # Split the WoO Map
        
    # the output will be the tiles of raster files with extreme value of salinity 
    # that limits the growth of the existing trees and kill the sapling directly
    
    ## Call the MesoFON and run with the newly created raster world
    # pause after 1 step (=3 months)
    # ~ MesoFON will pause and wait for the new raster files from Delft3D
    
    # receive the new trees location signal from MesoFON
    
    ## Rebuild the trees based on the tiles 

### End Loop

#%%
## Post Processing for the Delft3D Simulation
## Post Processing for the MesoFON Simulation  


