# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 14:44:53 2021

@author: sbe002
"""
import os
# import numpy as np
import numpy.ma as ma

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist

import sys
# print(sys.path)
sys.path.append('D:/Git/d3d_meso/FnD3D') # as this Func will be in the same folder, no longer needed
from d3d_meso_mangro import create_xyzwCellNumber, create_xyzwNodes

dir_testinput = os.path.join(r'D:/Git/d3d_meso/Model-Execute/D3DFM/FunnelMorphMF30')
file_nc_map = os.path.join(dir_testinput,'dflowfm','Delta_Schematized_funnel_net.nc')

# vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_map)
# ugrid_all = get_netdata(file_nc=file_nc_map)

# ugrid_all.mesh2d_face_x(mesh2d_face_nodes[:,2])

mesh_face_x = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_x')
mesh_face_x = ma.compressed(mesh_face_x)

mesh_face_y = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_y')
mesh_face_y = ma.compressed(mesh_face_y)

mesh_face_nodes = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_nodes')
# shp_mesh_face_nodes = mesh_face_nodes.shape

# mesh_face_nodes_compreessed = ma.getdata(mesh_face_nodes)

# a = mesh_face_x[mesh_face_nodes_compreessed[:,2]]


#%% Coba tes sama dengan bmi var yang mana
import numpy as np
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

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
## Initiate the BMI Delft3D
## Paths for dll-files and input-files for DFM
D3D_HOME = os.path.join(r'C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64')
MFON_HOME = os.path.join(r'Model-Execute/MesoFON')
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
# print ('DFM initialized')

## Tes korelasi cell number
# cell center coordinate: xz, yz
# node coordinate : xk, yk
model_dfm.get_var_shape('ba') # to get the size of the parameter

xk = model_dfm.get_var('xk') #Net node x coordinate {"shape": ["numk"]}
yk = model_dfm.get_var('yk')
#xz, yz, related with bl, s1, 
xz = model_dfm.get_var('xz') #waterlevel point / cell centre, x-coordinate (m) {"location": "face", "shape": ["ndx"]}
yz = model_dfm.get_var('yz') #y coordinate
#xzw, yzw related with cell number
xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
yzw = model_dfm.get_var('yzw') #y coordinate

model_dfm.get_var('ba')

#### calculate the cell number as in the position of xzw and yzw or ndxi
xyzw_cell_number = create_xyzwCellNumber(xzw,yzw,mesh_face_x,mesh_face_y)
xyzw_nodes = create_xyzwNodes(mesh_face_nodes,xyzw_cell_number)

#### Nodes coordinate sample script
# remember the cell number is not in sequence as in the grid,
# whereas it follows the xzw and yzw order
# Therefore, the rearranged nodes position based on the BMI-based cell number
# is given in the xyzw_nodes variable.

# The xz and yz are actually used by the DFM to store the bl, ba, s1, and etc.
# xz-yz and xzw-yzw are almost identic, except the addition boundary point in
# the end of the array.
# For example, if the length of xzw-yzw is 986 rows, and  xz-yz is 1018 rows.
# So that, in practical, rows 0-985 belong to the internal element and
# 986-1017 belong to the boundary cell.

# The script below provides the example on how to retrieve the nodes coordinate
# that is stored in xk and yk.

#### retrieve nodes coordinate in index[0] as indexed in xyzw_nodes
# If you want to check the cell number index the xyzw_cell_number in column 2.
position = 256
data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
print(xyzw_cell_number[position,2]+1) # add 1 to adjust the real cell number in netCDF
# access the nodes number and get the coordinate
print(xk[data-1], yk[data-1])# substracted to 1 in order to adjust the 0-based position in python




        

