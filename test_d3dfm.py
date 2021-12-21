# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 12:19:41 2021

@author: sbe002
"""

# =============================================================================
# Ini adalah script yang ditulis dari hasil VC dengan Uwe
# 
# Coba ditulis sesuai dengan model dia
# 
# Total script ada sekitar 800an lines, akan tetapi coba ditulis dulu semua aja.
# =============================================================================

# Call libraries needed for computation
from ctypes import *
import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import bmi.wrapper
import scipy.integrate as si
from scipy import integrate # kok dobel import??
import os
import pandas as pd
import glob
import sys
import faulthandler
faulthandler.enable()
import datetime
import pretty_errors # to make errors look readable, colorful, and easy to debug
pretty_errors.configure(
    separator_character = '*',
    filename_display    = pretty_errors.FILENAME_EXTENDED,
    line_number_first   = True,
    display_link        = True,
    lines_before        = 5,
    lines_after         = 2,
    line_color          = pretty_errors.RED + '> ' + pretty_errors.default_config.line_color,
    code_color          = '  ' + pretty_errors.default_config.line_color,
    truncate_code       = True,
    display_locals      = True
)
pretty_errors.blacklist('c:/python')

# specify paths of dll files and input files
D3D_HOME = r'E:\00.Delft3D_FM_exe\00 From Uwe\FM_Compiled_1.2.8.62485\x64'
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm', 'bin', 'dflowfm.dll')

workdir = r'E: ....'
config_file = os.path.join(workdir, 'dimr_config.xml')
mdu_file = os.path.join(workdir, 'dflowfm', '001.mdu')
figsavefolder = os.path.join(workdir, 'figures')

#Create target directory if none exists
# if not os.path.exists(figsavefolder):
#     os.mkdir(figsavefolder)
#     print("Directory", figsavefolder, "created")
# else:
#     print("Directory", figsavefolder, "already exists")

print("D3D_HOME    :" + D3D_HOME)
print("dimr_path   :" + dimr_path)
print("config_file :" + config_file)

#Add corrects locations to environment variable path
os.environment['PATH'] = os.path.join(D3D_HOME, 'share', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dflowfm', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dwaves', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'esmf', 'scripts') \
+ ";" + os.path.join(D3D_HOME, 'swan', 'scripts')

print("PATH" + os.environ['PATH'])

#Define DFM Wrapper
model_dfm = bmi.wrapper.BMIWrapper(engine=dflowfm_path, configfile=mdu_file)

#Define and initialise DIMR Wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()
print("Model Initialize")

#Retrieving important model variation from FlowFM
#source: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/dflowfm/packages/dflowfm_lib/include/bmi_get_var.inc
ndx = model_dfm.get_var('ndx') #number of flow nodes (internal boundary)
ndxi = model_dfm.get_var('ndxi') #number of internal flowcells (2D + 1D)
xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
yzw = model_dfm.get_var('yzw') #y coordinate of the center of gravity of the boxes
lnx = model_dfm.get_var('lnx') #total number of flow links between the boxes (internal boundary)
lnxi = model_dfm.get_var('lnxi') #number of links between within domain boxes (internal, 1D+2D)
ln = model_dfm.get_var('ln') #link matrix between adjacent boxes (node admin)
dx = model_dfm.get_var('dx') #distance between the centers of adjacent boxes (link length)
wu = model_dfm.get_var('wu') #width of the interface between adjacent boxes (link width)
ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area)
s0 = model_dfm.get_var('s0') #waterlevel    (m ) at start of timestep
s1 = model_dfm.get_var('s1') #waterlevel    (m ) at end   of timestep
hs = model_dfm.get_var('hs') #water depth at the end of timestep
bl = model_dfm.get_var('bl') #bottom level (m) (positive upward)

#Retrieving the setup on the vegetation and boundaries installation
#Properties here are stipulated for the trachytope 154 application computation
rnveg = model_dfm.get_var('rnveg') # vegetation density --> 3D plant density , 2D part is basis input (1/m2)
diaveg = model_dfm.get_var('diaveg') # vegetation diameter
Cdvegsp = model_dfm.get_var('Cdvegsp') # **this creates an issue** --> Uwe --> spatial plant Cdveg; maybe related to this paper: https://agupubs-onlinelibrary-wiley-com.tudelft.idm.oclc.org/doi/pdf/10.1002/2013JF002875
stemheight = model_dfm.get_var('stemheight')
botdep = model_dfm.get_var('bl') #bottom level (m) (positive upward) at the end of timestep --> redundant since I wrote it as in Uwe 
#botdep_adj = botdep[range(ndxi)]

#wetdry_thrshld = 0.05 # the minimum depth water depth in a cell to be considered wet (m)
#
#Add perturbations
#botdep = botdep+((np.random.rand(ndx)*0.04)-0.02)
#model_dfm.set_var('bl', botdep)
#
print("ndx" + str(ndx))
print("ndxi" + str(ndxi))
# print("yzw" + str(yzw))
print("lnxi" + str(lnxi))
print("ln" + str(ln))
# print("dx" + str(dx))
# print("wu" + str(wu))
# print("ba" + str(ba))
# print("s0" + str(s0))
# print("s1" + str(s1))
# print("hs" + str(hs))
# print("bl" + str(bl))
# print("rnveg" + str(rnveg))
# print("diaveg" + str(diaveg))
# print("Cdvegsp" + str(Cdvegsp))
# print("stemheight" + str(stemheight))
# print("botdep" + str(botdep))

#Boundaries for vegetation, the general area where it is possible to have vegetation growth
#
xbndmin = min(xzw[range(ndxi)]) #is this specified based on the cells where the water depth is below a certain level
xbndmax = max(xzw[range(ndxi)]) 
ybndmin = min(yzw[range(ndxi)])+8.0
ybndmax = max(yzw[range(ndxi)])-8.0

#Graphical outputs
make_veg_plot = 1
show_veg_plot = 1

#Triangulate face coordinates for plotting
face_triangulate = tri.Triangulation(xzw[range(ndxi)],yzw[range(ndxi)])
#Interpolate on the grid
#2D Plot
#define basic plotting routine for model output, (fld=variable, face_tiang, fsf=figsavefolder, i*timestep, lvls=, ttl=tit
def showfid(fld, face_triang, fsf,i,lvls,tl):
    plt.tricontourf(face_triang,fld)
    plt.title(ttl)
    plt.colorbar(shrink=0.3, extend = 'both', aspect=10)
    plt.axes().set_aspect('equal')
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    figsavepath = os.path.join(figsavefolder, ttl + '_tw' + str(i))
    plt.savefig(figsavepath, dpi=300, bbox_inches = 'tight')
    if show_veg_plot==1:
        plt.show()
    plt.close()

def showcross (face_triang, fld, fld1, fld2, fld3, i):
    x = np.arrange(0,ndx*10,10)
    plt.figure()
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.subplot(121)
    plt.plot(x, fld, 'tab:brown', label='Bed Level')
    plt.plot(x, fld1, 'c', label='Water Level')
    plt.plot(x, fld2, 'g--', label='Vegetation Height') #stem height
    plt.title('Bed Level and Water Level')
    plt.xlabel('Distance in the Cross Shore Direction')
    plt.ylabel('Elevation (m)')
    plt.grid(True)
    plt.legend(loc='bottom left')
    
    plt.subplot(122)
    plt.tricontourf(face_triang, fld3)
    plt.title('Diameter Layout')
    plt.colorbar(shrink=0.3, extend='both', aspect=10)
    plt.xlabel('Coordinates')
    plt.ylabel('Coordinates')
    figsavepath_cross=os.path.join(figsavefolder + 'Diameter Layout' + '_t-' + str(i))
    plt.savefig(figsavepath_cross, dpi=300, bbox_inches='tight')
    
    if show_veg_plot==1:
        plt.show()
    plt.close
    
#Set time parameters for coupled model
MF=100
tstep=(365*24*60*60)/MF #coupling yearly in seconds
# tstep=(90*24*60*60)/MF #coupling 3 monthly in seconds
vege_years = 20 #hydrodynamics #number of vegetation time steps
nt = 20
# nt = int(nt)

dt = 1.0 # time resolution (years)
mtpervt = (tstep*MF*vege_years)/(MF*nt) #seconds of model time per vegetation step
ncells = len(rnveg)
print('# of vegetation timesteps' + str(nt))
print('Vegetation timestep length' + str(mtpervt) + 'seconds')

#Mangrove Growth parameter
Seed = 0.05 #chance of establishment of seedlings in a grid cell (per yr)
# Seed = 0.05/(365*24/MF) #chance of establishment of seedlings in a grid cell (per hr MF)
Dmax = 140  #maximum diameter (cm)
P0 = 30 #initial mangrove density, max sapling recruits per year per plot, trees/100m2
D0 = 0.0137 #initial diameter, (1.37m @ breast height)
Max_age = 300 #Years
Hmax = 3500 #maximum height possible for the mangrove tree, cm
GG = 162 #species-specific growth constant, cm/year
# GG = 152.17/(36524/MF) #species-specific growth parameter, cm/year---cm/hr.MF
b2 = 48.04 #species-specific growth parameter, unitless
b3 = 0.172 #species-specific growth parameter, cm^(-1)
# amp = 2+0.15 #tidal amplitude + hlf of the wave height, m
amp = max(s0) #max water level
diam_def = 1.4 #default diameter stems, max, m
height_def = 35 #default height stems, max, m
nap_def = 30
Dif_veg = 0.2

#Inundation & Competition Stress Coefficient
aa = 4 #constant
bb = -8 #constant
cc = 0.5 #constant
d_const = -0.005 #constant in relation to the competition stress
#B_half = 835 #the value for the biomass, B when competition stress = 0.5 ~1044kg/m2

#Pneumatophore & Drag Coefficient Calculation
f = 30 #constant
Cdno = 0.005 #the value for the drag coefficient without the pressence of vegetation
e = 2 #the dimensional constant

####### Not complete

