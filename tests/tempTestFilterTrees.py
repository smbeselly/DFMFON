# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 17:27:42 2021

@author: sbe002
"""
#%% Import the libraries
import pandas as pd
import os
import numpy.ma as ma

#%% First initialization
#### NOTES: run the mesh reading first and then run BMI.
# Because it accesses the same dll for netCDF4.
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
from scipy import integrate
import faulthandler
faulthandler.enable()
import sys
# print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
# from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes   

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

#### Do the calculation
## Paths for dll-files and input-files for DFM
D3D_HOME = os.path.join(r'C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64')
MFON_HOME = os.path.join(r'D:/Git/d3d_meso/Model-Execute/MesoFON')
# workdir = os.path.join(r'D:\Git\d3d_meso\Model-Execute\D3DFM') # model rectangular
workdir = os.path.join(r'Model-Execute/D3DFM/FunnelMorphMF30') # model funnel with morpho

# Coba dari instalasi
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm','bin')
dflowfm_engine = os.path.join(dflowfm_path, 'dflowfm.dll')
config_file = os.path.join(workdir, 'FunnelMorphMF30.xml') # funnel morpho

mdu_file = os.path.join(workdir, 'dflowfm', 'FlowFM.mdu')
grid_file = os.path.join(workdir, 'dflowfm', 'Delta_Schematized_funnel_net.nc')
figsavefolder= os.path.join(r'Model-Out/D3DFM/Figs')

# Point to the data
dir_testinput = os.path.join(r'D:/Git/d3d_meso/Model-Execute/D3DFM/FunnelMorphMF30')
file_nc_map = os.path.join(dir_testinput,'dflowfm','Delta_Schematized_funnel_net.nc')

### Access the information from mesh
mesh_face_x = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_x')
mesh_face_x = ma.compressed(mesh_face_x)
mesh_face_y = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_y')
mesh_face_y = ma.compressed(mesh_face_y)
mesh_face_nodes = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_nodes')

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

#%% Let's get the job done
## TODO
# 3. iterasi untuk semua pohon dengan cara:
#   a. Filter pohon berdasarkan xk,yk
#   b. hitung species-dependent Htrunk, Dtrunk, Hpneu, Dpneu
#   c. Hitung drag coefficient representative
# 4. Lakukan pada masing2 cell yg mana cell sudah diseleksi berdasarkan WoO 

# We already have the function to point out the coordinates of the nodes
# Use that information to do the selection

## The trees location is in text file either as predefined or model results
# it as txt (csv) file with this following header:
# 'run','dbh_cm','Height_cm','Age','tick','Species','GeoRefPosX','GeoRefPosY'

# use pandas to access the file
MFON_HOME = os.path.join(r'D:/Git/d3d_meso/Model-Execute/MesoFON')
file_loc = os.path.join(MFON_HOME, 'instance_1','MF_Trees_Out.2021.Dec.08.14_26_54.txt')

read_data = pd.read_csv(file_loc)

## How to:
# 1 baca xyzw_cell_number, tunjukkan cell number di column [2]
# 2 access xk,yk
# 3 filter pohon berdasarkan xmin,xmax dan ymin,ymax

row = 5
position = xyzw_cell_number[row,2].astype(int)
nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
# Find the min max of each x,y coordinate
# create the list of x_min-x_max and y_min-y_max
x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]

## untuk test saja dummy value
x_range = [707150, 707200]
y_range = [9163130, 9163180]

# subsetting pandas 
#source: https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
# and https://www.geeksforgeeks.org/selecting-rows-in-pandas-dataframe-based-on-conditions/
read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                             (read_data['GeoRefPosX'] <= x_range[1])]
read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                    (read_data_subset['GeoRefPosY'] <= y_range[1])]
