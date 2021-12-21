# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Script to Coupling Delft3D and MesoFON

#%%
# Import the necessary packages
# import pywintypes
# import win32api
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
# see note: https://stackoverflow.com/questions/20554074/sklearn-omp-error-15-initializing-libiomp5md-dll-but-found-mk2iomp5md-dll-a

#%%
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

# engine and config setting MesoFON
# TODO
# mfon_path = os.path.join(MFON_HOME, 'src', 'meso_FON', 'DFM_Wrapper.java') ## maybe not necessary anymore
# mfon_config_file = os.path.join(MFON_HOME, 'batch', 'batch_params.xml')

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

## Define and initialise DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()
# print ('DFM initialized')

# ## Define MesoFON Wrapper
# model_mfon = bmi.wrapper.BMIWrapper(engine=mfon_path, config_file=mfon_config_file)
# model_mfon.initialize()
# print('MesoFON initialized')

#Retrieving important model variation from FlowFM
#source: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/dflowfm/packages/dflowfm_lib/include/bmi_get_var.inc
ndx = model_dfm.get_var('ndx') #number of flow nodes (internal boundary)
ndxi = model_dfm.get_var('ndxi') #number of internal flowcells (2D + 1D)
xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
yzw = model_dfm.get_var('yzw') #y coordinate
lnx = model_dfm.get_var('lnx') #total number of flow links between the boxes (internal boundary)
lnxi = model_dfm.get_var('lnxi') #number of links between within domain boxes (internal, 1D+2D)
ln = model_dfm.get_var('ln') #link matrix between adjacent boxes (node admin) 1D link (2,*) node   administration, 1=nd1,  2=nd2   linker en rechter celnr {"shape": [2, "lnkx"]}
dx = model_dfm.get_var('dx') #distance between the centers of adjacent boxes (link length)
wu = model_dfm.get_var('wu') #width of the interface between adjacent boxes (link width)
ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
s0 = model_dfm.get_var('s0') #waterlevel    (m ) at start of timestep {"location": "face", "shape": ["ndx"]}
s1 = model_dfm.get_var('s1') #waterlevel    (m ) at end   of timestep {"location": "face", "shape": ["ndx"]}
hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
bl = model_dfm.get_var('bl') #bottom level (m) (positive upward) {"location": "face", "shape": ["ndx"]}
xz = model_dfm.get_var('xz') # waterlevel point / cell centre, x-coordinate (m) {"location": "face", "shape": ["ndx"]}
yz = model_dfm.get_var('yz') # waterlevel point / cell centre, y-coordinate (m)
lncn = model_dfm.get_var('lncn')
LLkkk = model_dfm.get_var('LLkkk')


model_dfm.get_var_type('bl')
model_dfm.get_var_shape('xz') # to get the size of the parameter

#it seems like this is the coordinate of the nodes
xk = model_dfm.get_var('xk') #Net node x coordinate {"shape": ["numk"]}
yk = model_dfm.get_var('yk')

numk = model_dfm.get_var('numk')
is_dtint=model_dfm.get_var('is_dtint') # total time interval since last statistics reset. {"rank": 0}
is_maxvalsnd=model_dfm.get_var('is_maxvalsnd')	# [-] Integral values on flow nodes
is_sumvalsnd=model_dfm.get_var('is_sumvalsnd') # Integral values on flow nodes. {"location": "face", "shape": ["is_numndvals", "ndx"]}


# =============================================================================
# to check variables available in DFM
# for i in range(model_dfm.get_var_count()):
#     print(model_dfm.get_var_name(i))
# =============================================================================

# check time stepping
model_dfm.get_start_time(), model_dfm.get_current_time(), \
    model_dfm.get_end_time(), model_dfm.get_time_step()

#%% Test initialization by plotting the data 
# see link: https://github.com/openearth/notebooks/blob/master/delft3d%20bmi%20vrui.ipynb

# os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
# # see note: https://stackoverflow.com/questions/20554074/sklearn-omp-error-15-initializing-libiomp5md-dll-but-found-mk2iomp5md-dll-a

data = {
    'bl': model_dfm.get_var('bl'),
    's1': model_dfm.get_var('s1'),
    's0': model_dfm.get_var('s1').copy(),
    'xz': model_dfm.get_var('xz'),
    'yz': model_dfm.get_var('yz')
}

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
N = matplotlib.colors.Normalize(data['bl'].min(), data['bl'].max())
cmap = cmocean.cm.deep_r
axes[0].scatter(data['xz'], data['yz'], c=cmap(N(data['bl'])), 
                edgecolor='none')
# axes[0].scatter(data['xzw'], data['yzw'], c=data['bl'], cmap=cmap, edgecolor='none')
N = matplotlib.colors.Normalize(data['s1'].min(), data['s1'].max())
cmap = cmocean.cm.delta
axes[1].scatter(data['xz'], data['yz'], c=cmap(N(data['s1'])), 
                edgecolor='none')
# axes[1].scatter(data['xzw'], data['yzw'], edgecolor='none')


face_triangulate = tri.Triangulation(xzw[range(ndx)],yzw[range(ndx)])
plt.tricontourf(face_triangulate,s1)
#%%

# Start Run
model_dimr.update(3600)

realtime = 0

# Extract variables from DFM through BMI
is_dtint=model_dfm.get_var('is_dtint')
realtime = realtime+is_dtint
bedlvl=model_dfm.get_var('bl')
is_maxvalsnd=model_dfm.get_var('is_maxvalsnd')
is_sumvalsnd=model_dfm.get_var('is_sumvalsnd')

tau_mn=is_sumvalsnd[range(ndxi),0]/is_dtint


print('real time is ' + str(realtime) + ' s, or ')


model_dimr.finalize()

#%% coba seleksi mangrove berdasarkan cell
# =============================================================================
# 1. point out koordinat cell center (xz,yz)
# 2. berdasarkan (1) cari titik nodes terdekat (xk,yk)
# 3. iterasi untuk semua pohon dengan cara:
#   a. Filter pohon berdasarkan xk,yk
#   b. hitung species-dependent Htrunk, Dtrunk, Hpneu, Dpneu
#   c. Hitung drag coefficient representative
# 4. Lakukan pada masing2 cell yg mana cell sudah diseleksi berdasarkan WoO 

#source: https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
from scipy import spatial
ndx = model_dfm.get_var('ndx')
# A = np.zeros(shape=(ndx,2))
# A[:,0] = xz
# A[:,1] = yz

A = np.block([[xk],[yk]]).T #xk,yk is node points

pt = [xz[0], yz[0]]  # <-- the point to find xz,yz is the cell center points

# B = A[spatial.KDTree(A).query(pt)[1]] # <-- the nearest point    
index = spatial.KDTree(A).query(pt)
B = A[index[1]]

C = np.zeros((2,2))
C[0,0] = B[0]
C[1,0] = B[1]

# jika pt[0] > B[0] maka cari x dan yg lbh besar dari pt[0]
if pt[0]>B[0] and pt[1]>B[1] :
    C[0,1]=xz[1]
    C[1,1]=yz[1]
    
# =============================================================================
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
#Finalize the running
model_dimr.finalize()
print('model finished')
#%%
## Post Processing for the Delft3D Simulation
## Post Processing for the MesoFON Simulation  


