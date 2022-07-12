# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 19:05:21 2022

@author: sbe002

TODO belum nambah script untuk tambahan sediment karena biomass (daun jatuh, dll)
"""
# =============================================================================
# import datetime
# x1 = datetime.datetime.now()
# =============================================================================
#%% Input Folders and Files
# PROJ_HOME = r'D:\Git\d3d_meso'
# SYS_APP = r'D:\Git\d3d_meso/FnD3D'
# D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\2sebrian_20220518' #sal_veg-OK
# gdal_loc = r'C:\Users\sbe002\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
# JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'
# # Mangr_SHP = 'geserMangroveAgeMerged.shp'
# Mangr_SHP = 'Tip_Saplings_geser_2.shp'

# D3D_Model = 'FunnelMorphMF30_Adjusted_Saline_geser_2'
# D3D_Domain = 'Grid_Funnel_1_net.nc'
# MFON_Folder = 'MesoFON_20220506_noSurv'

PROJ_HOME = r'C:\Users\sbe002\Downloads\Research_Run\d3d_meso_run7_scenario_A_test'
SYS_APP = r'D:\Git\d3d_meso/FnD3D'
D3D_HOME = r'D:\Git\d3d_meso\Model-Execute\D3DFM\2sebrian_20220518' #sal_veg-OK
gdal_loc = r'C:\Users\sbe002\Anaconda3\envs\d3dfm_39\Lib\site-packages\osgeo_utils'
JAVA_Exe = r'C:\Users\sbe002\RepastSimphony-2.8\eclipse\jdk11\bin\java.exe'
# Mangr_SHP = 'Tip_Saplings_geser_2.shp' # tidak perlu karena akan randomly generated

D3D_Model = 'd3d_meso_run6_scenario_A'
D3D_Domain = 'Grid_Funnel_1_net.nc'
config_xml = 'd3d_meso_run6.xml'
mdu_file = 'FlowFM.mdu'
MFON_Folder = 'MesoFON_20220506_noSurv'

## Check the complete_model.jar file and change this source file
Sal_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_sal_'
Surv_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_surv_'
Excel_Source = r'C:\Users\brian\git\macro_FON_220111\meso_FON\tile_20_20_trees_input.xls'

#%% Import the necessary packages, set the file path, and input files

import os
import numpy as np
import numpy.ma as ma
import bmi.wrapper

import sys
print(sys.path)
sys.path.append(SYS_APP) # as this Func will be in the same folder, no longer needed
from d3d_prep_raster import d3dConcaveHull, d3dPolySHP, d3dCSV2ClippedRaster, d3dRaster2Tiles
# from d3d_meso_mangro import csv2ClippedRaster
from dfm_tools.get_nc import get_netdata, get_ncmodeldata
from dfm_tools.io.mdu import read_deltares_ini
from d3d_meso_mangro import calcAgeCoupling0, create_xzyzCellNumber, d3dNewRaster2Tiles
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
from d3d_mangro_seeds import seedling_prob, elim_seeds_surv, range_seed
import matplotlib.pyplot as plt
plt.close('all')

import glob
# from osgeo import ogr
# Supress/hide the warning
np.seterr(invalid='ignore')

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d
from xlrd import open_workbook
from xlutils.copy import copy
import shutil
from pathlib import Path
from osgeo import gdal, gdalconst
import zipfile
import send2trash
import sys
import fileinput
import re
from dateutil import parser
import datetime
import copy

## Set the paths for dll-files and input-files for DFM
PROJ_HOME = os.path.join(PROJ_HOME)
D3D_HOME = os.path.join(D3D_HOME)
MFON_HOME = os.path.join(PROJ_HOME,'Model-Execute','MesoFON')
D3D_workdir = os.path.join(PROJ_HOME,'Model-Execute','D3DFM',D3D_Model)
MFON_JAR = os.path.join(MFON_HOME, MFON_Folder,'complete_model.jar')
MFON_LocalBatchRunner = os.path.join(MFON_HOME,'local_batch_run.properties')
gdal_path = os.path.join(gdal_loc)
JAVAREP = os.path.join(JAVA_Exe)

MFON_Exchange = os.path.join(PROJ_HOME,'Model-Exchange')
if not os.path.exists(MFON_Exchange):
    os.makedirs(MFON_Exchange)
MFON_Env = os.path.join(MFON_Exchange, 'MesoFON-Env')
if not os.path.exists(MFON_Env):
    os.makedirs(MFON_Env)
MFON_Trees = os.path.join(MFON_Exchange, 'MesoFON-Trees')
if not os.path.exists(MFON_Trees):
    os.makedirs(MFON_Trees)
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')
if not os.path.exists(MFON_OUT):
    os.makedirs(MFON_OUT)
DFM_OUT = os.path.join(PROJ_HOME,'Model-Out','D3DFM')
if not os.path.exists(DFM_OUT):
    os.makedirs(DFM_OUT)
dir_out = os.path.join(MFON_Exchange, 'Initialization')
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
figsavefolder= os.path.join(PROJ_HOME,'Model-Out','Figures')
if not os.path.exists(figsavefolder):
    os.makedirs(figsavefolder)
seedlings_out = os.path.join(PROJ_HOME,'Model-Out','MesoFON', 'Seedlings')
if not os.path.exists(seedlings_out):
    os.makedirs(seedlings_out)
botdepth_out = os.path.join(DFM_OUT, 'Bottom_Depth')
if not os.path.exists(botdepth_out):
    os.makedirs(botdepth_out)

## settings

EPSG_Project = 32749 # EPSG code for WGS84/ UTM Zone 49S (Porong case study)
netcdf_domain = os.path.join(D3D_workdir, 'dflowfm', D3D_Domain)
x_res = 10
y_res = 10
no_data_val = -999.0
tile_size_x = 200 # it is in meter 
tile_size_y = 200 # which reads info in pixel
species_name = 'Avicennia_marina'
LLWL = -1.2
coupling_period = 90 #days, actually 90 days
MorFac = 30 
woo_inun = 3 # inundation free period (days)
no_data_val = -999.0
limit_seed = [706574,707031,9163156,9163511] # [xmin, xmax, ymin, ymax] the area allowed for seedlings development

#%% Import the nc file, create hull, and build poly
nc_in = netcdf_domain
LLWL_CH = LLWL - 0.5 # to make correction for concavehull process

ugrid_all = get_netdata(file_nc=nc_in)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# For instance the LLWL value is - 1.5 + 0.5m
matrix_llwl = np.where(matrix[:,2] < LLWL_CH , -999, matrix[:,2])
matrix = np.column_stack([matrix[:,:2],matrix_llwl])
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])
# create matrix x,y for hull
mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T 

#%% From Poly create SHP

projection = EPSG_Project
dir_out = dir_out
out_poly_name = 'CH_'

#%% Create raster from bathimetri .csv with gdal_grid
# Create .csv file
concave_path = dir_out
concave_name = 'CH_bathy_'
EPSG_coord = EPSG_Project
x_res = x_res
y_res = y_res
no_data_val = no_data_val

shp_clip = out_poly_name # this assume that the shp file similar with one
                        # build from d3dPolySHP
                        # It will directly refer to the shapefile in the same folder
                        # with 'CH_.shp'
affix = '_clipped'
# create matrix, in this case we can create a custom matrix from data reading
# for instance, selection of data points of WoO
# don't forget to adjust the shapefile from the allowable WoO

ugrid_all = get_netdata(file_nc=nc_in)#,multipart=False)
matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
# clean data from -999 and nan value
matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
matrix = (matrix[~np.isnan(matrix).any(axis=1)])

d3dCSV2ClippedRaster(concave_path, concave_name, EPSG_coord, matrix, x_res, y_res, no_data_val, shp_clip, affix)

#return to home
os.chdir(PROJ_HOME)
#%% Tile the raster and filter + delete nodata tiled raster
out_path = dir_out
output_filename = "tile_"

ras_clip = os.path.join(out_path+str('\\')+concave_name+affix+'.tif')

tile_x = tile_size_x # it is in meter 
tile_y = tile_size_y # which reads info in pixel

d3dRaster2Tiles(out_path, output_filename, ras_clip, tile_x, tile_y,CreateSHP=True)

import gc
gc.collect() # to clear memory of variables in python after doing del(variables)

#%% Run the DFM until reach the condition where accommodation space is available
### Initiate the BMI
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
### Access the variables from BMI.get_var
xk = model_dfm.get_var('xk') #Net node x coordinate {"shape": ["numk"]}
yk = model_dfm.get_var('yk')
#xz, yz, related with bl, s1, 
xz = model_dfm.get_var('xz') #waterlevel point / cell centre, x-coordinate (m) {"location": "face", "shape": ["ndx"]}
yz = model_dfm.get_var('yz') #y coordinate
#xzw, yzw related with cell number
xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
yzw = model_dfm.get_var('yzw') #y coordinate

### calculate the cell number as in the position of xzw and yzw or ndxi
ugrid_all = get_netdata(file_nc=grid_file)
xzyz_cell_number = create_xzyzCellNumber(xz, yz, model_dfm, ugrid_all)

#%% routine to run the DFM and pause when reach January to check the accommodation space

### get the reference date and starttime of model
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

### Loop the Coupling
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

#%% Run the DFM 

# initiate an empty array
bed_level = np.empty((len(xz),0)) 
water_level = np.empty((len(xz),0)) 
coupling_ntime_is = 0
while bed_level.size == 0:
    # current simulation time is
    cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) 
    # current simulation time multiply by MF
    cursecMF = cursec*MorFac
    # current month based on simulation time (with MF)
    curyr = ((refdatet+cursecMF).strftime('%Y'))
    curmonth = ((refdatet+cursecMF).strftime('%m'))
    timeisnow = datetime.datetime.now()   
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        model_dimr.update()       
        ## check to apply seedling establishment
        if curmonth == '01':
            print('prepare for natural seedlings establishment')
            bl = model_dfm.get_var('bl') # bottom level
            bed_level = np.append(bed_level, np.reshape(bl,(len(bl),1)), axis=1)
            ## get the water level value for the WoO
            s1 = model_dfm.get_var('s1') # water level
            # store the maximum water level per time step in column wise
            water_level = np.append(water_level, np.reshape(s1,(len(s1),1)), axis=1)
                 
        else:
            print(curyr, '/', curmonth, 'no seedlings establishment')
        
        dts = model_dfm.get_time_step()
        t=t+dts
        
        print('Coupling ',str(coupling_ntime_is), 'run ', t, '/', coupling_period_model)
    
    np.savetxt(os.path.join(botdepth_out, 'Coupling_'+str(coupling_ntime_is)+'.txt'), 
               model_dfm.get_var('bl'), delimiter=",")
    
    coupling_ntime_is += 1
    # save bottom depth as text to facilitate the ongoing simulation visualisation 
    
    try:
        # get the filtered bed_level
        bed_level_check = range_seed(bed_level, limit_seed, LLWL, no_data_val, xz, yz)
        
        if bed_level_check.size != 0:
            # init_seed_path = os.path.join(MFON_Exchange, 'init seed coupling'+str(ntime+1))
            # if not os.path.exists(init_seed_path):
            #     os.makedirs(init_seed_path)
            # surv_val_name = 'init_accm_space'
            # csv2ClippedRaster(init_seed_path, bed_level_check, surv_val_name, 
            #                   x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
            # #return to home
            # os.chdir(PROJ_HOME)
            from shapely.geometry import Polygon, Point, shape
            from scipy.spatial import ConvexHull, convex_hull_plot_2d
            hulls = ConvexHull(bed_level_check[:,:2])
            polylist = []
            for idx in hulls.vertices: #Indices of points forming the vertices of the convex hull.
                polylist.append(bed_level_check[:,:2][idx])
            poly_seed=Polygon(polylist)
            
            # poly_seed = Polygon(list(zip(bed_level_check[:,0], bed_level_check[:,1])))
            # poly_seed = Polygon(list(zip(hulls.simplices[:,0], hulls.simplices[:,1])))
        
        else:
            bed_level = np.empty((len(xz),0)) 
            water_level = np.empty((len(xz),0)) 
    except:
        # another confirmation that no seedlings establishment
        print(curyr, '/', curmonth, 'no seedlings establishment', 'Coupling', coupling_ntime_is)

#%% Calculate the seedlings establishment        
## The calculation or determination of the seedling establishment has ended
#  now to calculate the random distribution of the seedlings based on the
#  accommodation space
## define the randomizer function with numpy
# =============================================================================
# ## plot the poly_seed for checking
# 
# plt.plot(bed_level_check[:,0], bed_level_check[:,1], 'o')
# for simplex in hulls.simplices:
#     plt.plot(bed_level_check[:,:2][simplex, 0], bed_level_check[:,:2][simplex, 1], 'k-')
  
# p = gpd.GeoSeries(poly_seed)
# p.plot()

# check the area to shp
# gdf = gpd.GeoSeries([poly_seed])
# gdf.to_file("D:/test.shp")
# =============================================================================

import random
def poly_rand_pts(poly, num_points):
    min_x, min_y, max_x, max_y = poly.bounds
    points = []
    
    while len(points) < num_points:
        random_point = Point([random.uniform(min_x, max_x),
                              random.uniform(min_y, max_y)])
        if (random_point.within(poly)):
            points.append(random_point)
    
    return points

## get the area of the polygon first
# from pyproj import CRS
# def get_poly_area(poly, EPSG_Project):
#     # extract the outer ring of the polygon
#     import pyproj
#     w,s,e,n = poly.bounds
#     # pa = pyproj.Proj(init='EPSG:'+str(EPSG_Project), lat_1 = s, lat_2 = n)
#     pa = pyproj.Proj(CRS.from_epsg(EPSG_Project), lat_1 = s, lat_2 = n)
#     x,y = poly.exterior.coords.xy
#     # projection
#     x, y = pa(list(x[:-1]), list(y[:-1]))
#     #create a new geojson of the new projection polygon
#     cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
#     #geojson -> shapely polygon, 
#     #then use the polygon's predefined function to calculate area.
#     area = shape(cop).area # in m^2
    
#     return area

## define random seedlings position
# find number of seedlings per area: 0.03/m^2 as in Porong
# area_seed = get_poly_area(poly_seed, EPSG_Project)
area_seed = poly_seed.area
seeds_are = round(area_seed*0.03) 
seeds_pt = poly_rand_pts(poly_seed,seeds_are)

xv = np.array([point.x for point in seeds_pt])
yv = np.array([point.y for point in seeds_pt])

# =============================================================================
# # plot of polygon area and random points
# p = gpd.GeoSeries(poly_seed)
# p.plot()
# plt.scatter(xv,yv, c='red')
# =============================================================================

### Calculate the WoO value
from d3d_meso_mangro import Calc_WoO
med_h_wl, surv_val = Calc_WoO(water_level, model_dfm, MorFac, coupling_period, woo_inun, LLWL)

## Get the filter of the seedlings based on the surv-val

seeds_pd = pd.DataFrame(data=np.column_stack((xv,yv)), columns=['GeoRefPosX', 'GeoRefPosY'])
seeds_pd['height'] = 1
seeds_pd['Age'] = 2
seeds_pd_2 = seedling_prob(seeds_pd, xzyz_cell_number, 
                                    ugrid_all, surv_val)
seeds_pd_filt = elim_seeds_surv(seeds_pd_2, xzyz_cell_number, ugrid_all, surv_val)

#%% Run again for the next two years before the seedlings can be embedded 
# in the DFM model
coupling_ntime_after_seeds = 7 # it represents 7 couplings or 1 year and 9 months
# this run without water level and salinity collection

for ntime in range(int(coupling_ntime_after_seeds)):
    # do the calculation for each coupling_ntime
    print('Start the coupling',str(ntime+coupling_ntime_is),'computation')
    cursec = datetime.timedelta(seconds=model_dfm.get_current_time())
    # current simulation time multiply by MF
    cursecMF = cursec*MorFac
    # current month based on simulation time (with MF)
    curyr = ((refdatet+cursecMF).strftime('%Y'))
    curmonth = ((refdatet+cursecMF).strftime('%m'))
    print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS')))
    
    timeisnow = datetime.datetime.now()   
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        model_dimr.update()
        dts = model_dfm.get_time_step()
        t=t+dts
        
        print('Coupling ',str(ntime+coupling_ntime_is), 'run ', t, '/', coupling_period_model)

        cursec = datetime.timedelta(seconds=model_dfm.get_current_time())
        # current simulation time multiply by MF
        cursecMF = cursec*MorFac
        # current month based on simulation time (with MF)
        curyr = ((refdatet+cursecMF).strftime('%Y'))
        curmonth = ((refdatet+cursecMF).strftime('%m'))
        print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
    
    np.savetxt(os.path.join(botdepth_out, 'Coupling_'+str(ntime+coupling_ntime_is)+'.txt'), 
               model_dfm.get_var('bl'), delimiter=",")    
    timeisend = datetime.datetime.now()
    print('Runtime', 'Coupling ',str(ntime+coupling_ntime_is), 'is', timeisend-timeisnow)

current_coupling_time = ntime+coupling_ntime_is

#%% Create and prepare for MesoFON Coupling
# run DFM for 1 coupling step and collect env information
# current simulation time is
cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) 
# current simulation time multiply by MF
cursecMF = cursec*MorFac
# current month based on simulation time (with MF)
curyr = ((refdatet+cursecMF).strftime('%Y'))
curmonth = ((refdatet+cursecMF).strftime('%m'))
print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
#initiating empty array for residual current calculation
res_x = np.empty((len(xz),0))
res_y = np.empty((len(xz),0))   
water_level = np.empty((len(xz),0)) 
salinity = np.empty((len(xz),0))
ntime = current_coupling_time

timeisnow = datetime.datetime.now()   
t=0 # since the time step in DFM is flexible, therefore use this approach.
while t<coupling_period_model:
    model_dimr.update()
    s1 = model_dfm.get_var('s1') # water level
    sa1 = model_dfm.get_var('sa1') # salinity
    # store the maximum water level per time step in column wise
    water_level = np.append(water_level, np.reshape(s1,(len(s1),1)), axis=1)
    salinity = np.append(salinity, np.reshape(sa1,(len(sa1),1)), axis=1)
    
    ## check to apply seedling establishment
    if curmonth == '01' and curyr_check == 0:
        print('prepare for seedlings establishment')
        res_x, res_y = collect_res( model_dfm, res_x, res_y)
             
    elif curyr_check != 0:
        if curmonth == '01' and curyr != curyr_check:
            print('prepare for seedlings establishment')
            res_x, res_y = collect_res( model_dfm, res_x, res_y)
        else:
            print(curyr, '/', curmonth, 'no seedlings establishment')
    
    dts = model_dfm.get_time_step()
    t=t+dts
      
    print('Coupling ',str(ntime+1), 'run ', t, '/', coupling_period_model)

    cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) 
    # current simulation time multiply by MF
    cursecMF = cursec*MorFac
    # current month based on simulation time (with MF)
    curyr = ((refdatet+cursecMF).strftime('%Y'))
    curmonth = ((refdatet+cursecMF).strftime('%m'))
    print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
  
timeisend = datetime.datetime.now()
print('Runtime', 'Coupling ',str(ntime+1), 'is', timeisend-timeisnow)
curyr_check = curyr

if res_x.size != 0 :
    res_of_x = calculate_residual(res_x, model_dfm.get_time_step(), xz).reshape((len(res_x),1))
    res_of_y = calculate_residual(res_y, model_dfm.get_time_step(), xz).reshape((len(res_y),1))
    # create matrix (array)
    residual_is = np.hstack((res_of_x,res_of_y))
    print('calculate residual current in coupling',str(ntime+1))
else:
    residual_is = []
med_sal = np.median(salinity, axis=1)
print('Calculate median salinity in coupling',str(ntime+1))
#%% After DFM run processing
### Calculate the WoO value
read_data = seeds_pd_filt
age_coupling = read_data['Age']

# For loop for all of the cell number
index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
# initiate calculate subsurface contribution of mangrove roots
# bl_val = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
# addition_bl = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_val[-1]
# # model_dfm.get_var('bl')[bl_val.shape[0]-1:-1]
# bl_val = np.append(bl_val, addition_bl)


# from d3d_meso_mangro import Calc_WoO
med_h_wl, surv_val = Calc_WoO(water_level, model_dfm, MorFac, coupling_period, woo_inun, LLWL)

### 4. Convert from data point to raster environment 
# 4.1 Create raster from the surv-val
surv_val_raster = np.column_stack((xz,yz,surv_val))
concave_path = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
if not os.path.exists(concave_path):
    os.makedirs(concave_path)
surv_val_name = 'CH_surv_val'
csv2ClippedRaster(concave_path, surv_val_raster, surv_val_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
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

data_4shp = copy.copy(read_data)
data_4shp['Height_cm'] = data_4shp['height']*100
data_4shp['dbh_cm'] = 3
createPointSHP(data_4shp, age_coupling, concave_path, EPSG_Project) 

# 6.2. Tile the Master Trees and Save as XLS Input File
shp_source = os.path.join(concave_path, 'coupling'+str(ntime+1)+'.shp') # location of the source tree and the master tree shp
folder_loc = dir_out # location of the master tiles
file_tile = os.path.join(dir_out,'tile_*.shp' )
save_tiled_trees = os.path.join(MFON_Trees,'coupling'+str(ntime+1)) # location to save the tiled shp
if not os.path.exists(save_tiled_trees):
    os.makedirs(save_tiled_trees)
  
# clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name, a0, b0, a137, b137)
New_clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name)
    
### 8. Create the env raster and tile raster files
save_tiled_env = os.path.join(MFON_Env, 'coupling'+str(ntime+1))
if not os.path.exists(save_tiled_env):
    os.makedirs(save_tiled_env)
gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
# calc_surv = '"0*(A==-999)"'
# calc_sal = '"0*(A>-999)+35"'
# val_no_data_surv = 0
# val_no_data_sal = 60

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

# 9. Preparation for the MesoFON 
# 9.1 Copy and Place the payload.jar, and Prepare the Tiled xls to the designated folder

file_tile_xls = os.path.join(save_tiled_trees,'tile_*_input.xls')

# Delete the existing the execution folder tile_* to prepare for the new one
for filepate in glob.iglob(os.path.join(MFON_HOME,'tile_*')):
    # shutil.rmtree(filepath) # delete recursive
    # print(filepath)
    try:
        send2trash.send2trash(filepate)
    except OSError as e:
        print("Error: %s : %s" % (filepate, e.strerror))

# 9.2 Create new folders as the existing xls files and extract MesoFON jar file
for filepath in glob.iglob(file_tile_xls): # looping for all with trees affix
    # print(filepath)
   xls_folder = Path(filepath).stem
   os.makedirs(os.path.join(MFON_HOME,xls_folder))
   with zipfile.ZipFile(MFON_JAR, 'r') as zip_ref:
       zip_ref.extractall(os.path.join(MFON_HOME,xls_folder))
   shutil.copy(MFON_LocalBatchRunner, os.path.join(MFON_HOME,xls_folder))

# 9.3 Replace the content in Unrolled_Param, batch_params.xml, and parameters.xml
# function to replace the content of the file
# source = https://www.delftstack.com/howto/python/python-replace-line-in-file/
def replacement(file, previousw, nextw):
   for line in fileinput.input(file, inplace=1):
       line = line.replace(previousw, nextw)
       sys.stdout.write(line)
# The Looping to replace the variables in the params files
for filepatf in glob.iglob(os.path.join(MFON_Exchange,'Initialization','tile_*.tif')):
    Init_Rasters = save_tiled_env
    Init_Trees = save_tiled_trees
    # Create Update Path
    Sal_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_sal_').replace("\\",'/')
    Surv_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_surv_').replace("\\",'/')
    Excel_Update = os.path.join(Init_Trees,Path(filepatf).stem+'_trees_input.xls').replace("\\",'/')
    # Point the Files that need to be updated
    unrolledParam = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','unrolledParamFile.txt')
    batchParam = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','scenario.rs','batch_params.xml')
    param_xml = os.path.join(MFON_HOME,Path(filepatf).stem+'_trees_input','scenario.rs','parameters.xml')
    # Replace the content of the file
    # unrolledParamFile.txt
    replacement(unrolledParam, Surv_Source, Surv_Update)       
    replacement(unrolledParam, Sal_Source, Sal_Update)
    replacement(unrolledParam, Excel_Source, Excel_Update)
    # batch_params.xml
    replacement(batchParam, Surv_Source, Surv_Update)       
    replacement(batchParam, Sal_Source, Sal_Update)
    replacement(batchParam, Excel_Source, Excel_Update)
    # parameters.xml
    replacement(param_xml, Surv_Source, Surv_Update)       
    replacement(param_xml, Sal_Source, Sal_Update)
    replacement(param_xml, Excel_Source, Excel_Update)

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
    all_df.append(df)

Concat_table = pd.concat(all_df)

# 12.1. drop tick 0 year, only take 0.25
Concat_table = Concat_table[Concat_table.tick > 0]
run_is = 'Coupling_'+str(ntime+1) # change this with the real name

# 12.2. use spatial in scipy to match the x,y of the mangroves and the age information.
master_trees = gpd.read_file(os.path.join(concave_path+str('\\')+Path(concave_path).stem+'.shp'))
age_coupling = calcAgeCoupling0(Concat_table, master_trees) # prepare the age_coupling for the Master Trees of Each Coupling procedure
Concat_table['Age'] = age_coupling #update age after MesoFON run
Concat_table.drop(['tick'], inplace=True, axis=1)
Concat_table = Concat_table.reset_index(drop=True) # reset the index

# 12.3. Concatenated table is saved as txt file
Concat_table.to_csv(os.path.join(MFON_OUT_compile, run_is+'.txt'), sep=',', index=False, header=True)
# save bottom depth as text to facilitate the ongoing simulation visualisation 
np.savetxt(os.path.join(botdepth_out, run_is+'.txt'), model_dfm.get_var('bl'), delimiter=",")

### 13. Read the compile txt file to prepare for the next iteration
read_data = Concat_table  

# initiate calculate subsurface contribution of mangrove roots
bl_val = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
addition_bl = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_val[-1]
# model_dfm.get_var('bl')[bl_val.shape[0]-1:-1]
bl_val = np.append(bl_val, addition_bl)

# # calculate below ground biomass contribution to bed level
# bl_mangro = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
# addition_bl_mangro = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_mangro[-1]
# bl_mangro = np.append(bl_mangro, addition_bl_mangro)
# # calculate the delta after every coupling 
# bl_delta = bl_mangro - bl_val
# # save current calcLevelCell value for the next iteration
# bl_val = bl_mangro.copy()

timeislast = datetime.datetime.now()
print('End of coupling',str(ntime+1))
print('Date now with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 

#%% This part is to continue to the calculation of DFM with Cdveg and Mangrove
Mangr_SHP = concave_path+str('\\')+Path(concave_path).stem+'.shp'
#### Read Master Trees
master_trees = gpd.read_file(os.path.join(MFON_Trees, 'Master-Trees', Mangr_SHP))

#%% Read the compiled tree from the MesoFON Initialization Run and calculate drag_coefficient

MFON_OUT_compile = os.path.join(MFON_OUT,'Compile')

read_data = pd.read_csv(os.path.join(MFON_OUT_compile,'Coupling_'+str(ntime+1)+'.txt'))

# use spatial in scipy to match the x,y of the mangroves and the age information.
# age_coupling0 = calcAgeCoupling0(read_data, master_trees)
# age_coupling = calcAgeCoupling0(read_data, master_trees)
age_coupling = read_data['Age']

# For loop for all of the cell number
index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
drag_coeff = initCalcDraginLoopCdveg(xzyz_cell_number, model_dfm, ugrid_all, 
                                     index_veg_cel, read_data)
# assume that the boundary flow nodes are located at the end of array
addition = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + 0.005  
addition_veg = copy.copy(addition)*0
# drag_coeff = np.append(drag_coeff, addition)*ind
drag_coeff = np.append(drag_coeff, addition)

# initiate calculate subsurface contribution of mangrove roots
# bl_val = calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
# addition_bl = np.zeros((model_dfm.get_var('ndx')-model_dfm.get_var('ndxi'))) + bl_val[-1]
# # model_dfm.get_var('bl')[bl_val.shape[0]-1:-1]
# bl_val = np.append(bl_val, addition_bl)

#%% Loop the Coupling
# total couplings to be calculated
# coupling_ntime_after = coupling_ntimeUse - coupling_ntime_is - coupling_ntime_after_seeds
# coupling_ntime_run = coupling_ntimeUse - coupling_ntime_after_seeds -1
## Check the complete_model.jar file and change this source file
Sal_Source = 'coupling{}'.format(current_coupling_time+1)
Surv_Source = Sal_Source
Excel_Source = Sal_Source

coupling_ntime_run = coupling_ntimeUse - coupling_ntime_is - coupling_ntime_after_seeds
# current_coupling_time        
for ntime in range(int(coupling_ntime_run)):
    ntime += current_coupling_time+1
    # do the calculation for each coupling_ntime
    print('Start the coupling',str(ntime+1),'computation')
    ### 1.1. run the DFM all the simulation time within ntime
    # find cells that have vegetation
    index_veg_cel = index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data)
    # predefine the pandas dataframe of the mangrove positions
    list_read_subset = list_subsetCdveg(ugrid_all, xzyz_cell_number, index_veg_cel, read_data)
    
    rnveg_coeff, diaveg_coeff, stemheight_coeff = definePropVeg(xzyz_cell_number, 
                            model_dfm, ugrid_all, index_veg_cel, read_data, addition_veg)

    model_dfm.set_var('rnveg',rnveg_coeff)
    model_dfm.set_var('diaveg',diaveg_coeff)
    model_dfm.set_var('stemheight',stemheight_coeff)
    
    # update the variable with new value taken from previous drag calculation
    model_dfm.set_var('Cdvegsp',drag_coeff)
    water_level = np.empty((len(xz),0)) 
    salinity = np.empty((len(xz),0)) 
    #https://www.delftstack.com/howto/numpy/python-numpy-empty-array-append/
    
    # calculate mangroves' subsurface contribution
    try:
        bl_model = model_dfm.get_var('bl') + bl_delta # the correct script
        model_dfm.set_var('bl',bl_model)
    except:
        print('no bed level update, model initiation')
        
    try:
        list_seed2sapl.append(seedling_finalpos)
        seed2sapl = pd.concat(list_seed2sapl, axis=0)
        # add age to the appended pandas
        seed2sapl['Age'] = seed2sapl['Age']+add_seeds_age
        
        list_seed2sapl = []
        list_seed2sapl.append(seed2sapl)
        
        check_as_saplings = seed2sapl[seed2sapl['Age'] >= datetime.timedelta(days = 730)]

        if check_as_saplings.size > 0:
            print('Transfer seedlings as saplings at coupling', str(ntime+1))
            # select seedlings that have been transformed to saplings
            now_as_saplings = check_as_saplings.copy()
            now_as_saplings.loc[:, ['Age']] = (now_as_saplings['Age']/datetime.timedelta(days=365))
            now_as_saplings['dbh_cm'] = 3
            now_as_saplings['Height_cm'] = 100
            
            now_as_saplings.drop(['ParentPosX', 'ParentPosY','row',
                                  'Created Time Stamp', 'CrownSurfaceArea_m2'], 
                                 inplace=True, axis=1)
            
            # update the read data with new saplings as mature mangroves
            read_data = pd.concat([read_data, now_as_saplings], axis=0, ignore_index=True)
            
            # reset list to use the filterd list from current selection
            # filter for less than 730
            filter_seeds2apl = seed2sapl[seed2sapl['Age'] <= datetime.timedelta(days = 730)]
            # list_seed2sapl = []
            seed2sapl = filter_seeds2apl.copy()
            list_seed2sapl = []
            list_seed2sapl.append(seed2sapl)
            
            print('Total number of saplings:', now_as_saplings.shape[0])
            
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

    timeisnow = datetime.datetime.now()   
    t=0 # since the time step in DFM is flexible, therefore use this approach.
    while t<coupling_period_model:
        model_dimr.update()
        s1 = model_dfm.get_var('s1') # water level
        sa1 = model_dfm.get_var('sa1') # salinity
        # store the maximum water level per time step in column wise
        water_level = np.append(water_level, np.reshape(s1,(len(s1),1)), axis=1)
        salinity = np.append(salinity, np.reshape(sa1,(len(sa1),1)), axis=1)
        
        ## check to apply seedling establishment
        if curmonth == '01' and curyr_check == 0:
            print('prepare for seedlings establishment')
            res_x, res_y = collect_res( model_dfm, res_x, res_y)
                 
        elif curyr_check != 0:
            if curmonth == '01' and curyr != curyr_check:
                print('prepare for seedlings establishment')
                res_x, res_y = collect_res( model_dfm, res_x, res_y)
            else:
                print(curyr, '/', curmonth, 'no seedlings establishment')
        
        dts = model_dfm.get_time_step()
        t=t+dts
        
        if t % model_dfm.get_time_step() == 0:
            # calculate the drag coefficient
            drag_coeff = newCalcDraginLoopCdveg(model_dfm,xzyz_cell_number, 
                                                index_veg_cel,list_read_subset)
            drag_coeff = np.append(drag_coeff, addition)
            # update the variable with new value
            model_dfm.set_var('Cdvegsp',drag_coeff)
        
        print('Coupling ',str(ntime+1), 'run ', t, '/', coupling_period_model)

        cursec = datetime.timedelta(seconds=model_dfm.get_current_time()) #https://stackoverflow.com/questions/775049/how-do-i-convert-seconds-to-hours-minutes-and-seconds
        # current simulation time multiply by MF
        cursecMF = cursec*MorFac
        # current month based on simulation time (with MF)
        curyr = ((refdatet+cursecMF).strftime('%Y'))
        curmonth = ((refdatet+cursecMF).strftime('%m'))
        print('simulation date with MF is', ((refdatet+cursecMF).strftime('%Y%m%d %HH:%MM:%SS'))) 
      
    timeisend = datetime.datetime.now()
    print('Runtime', 'Coupling ',str(ntime+1), 'is', timeisend-timeisnow)
    curyr_check = curyr
    
    if res_x.size != 0 :
        res_of_x = calculate_residual(res_x, model_dfm.get_time_step(), xz).reshape((len(res_x),1))
        res_of_y = calculate_residual(res_y, model_dfm.get_time_step(), xz).reshape((len(res_y),1))
        # create matrix (array)
        residual_is = np.hstack((res_of_x,res_of_y))
        print('calculate residual current in coupling',str(ntime+1))
    else:
        residual_is = []
    med_sal = np.median(salinity, axis=1)
    # nozero = np.ma.masked_equal(salinity, 0)
    # med_sal = np.ma.median(nozero, axis=1)
    print('Calculate median salinity in coupling',str(ntime+1))
    
    # 1.3. Prepare pandas for list of the seedling_finalpos
    
    # check condition
    if len(residual_is) != 0:
        print('calculate seedlings position for year', curyr )
        seedling_finalpos = seedling_dispersal(xzyz_cell_number, index_veg_cel, ugrid_all, read_data, 
                               med_sal, residual_is, model_dfm, reftime, cursec)
        # check seedlngs duplicate
        bool_series = seedling_finalpos.duplicated(keep='first')
        seedling_finalpos = seedling_finalpos[~bool_series]
        
        seedling_finalpos.to_csv(os.path.join(seedlings_out, 
                                'Seedling_Coupling_'+str(ntime+1)+'.txt'), 
                                 sep=',', index=False, header=True)

    else:
        print(curyr, '/', curmonth, 'no seedlings establishment')
    
    
    ### 2. convert the water_level from each time_step to each day
    # calculate the x axis of the array of the default run
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
    # wl_shape = water_level.shape
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
    
    ### 3. Calculate the probability based on the WoO and store the value in each cell
    # find median value of the h_wl for each cell number
    med_h_wl = np.median(h_wl, axis=1) 
    # find median value of salinity for each cell number
    med_sal = med_sal
    
    # calculate WoO probability value from h_wl (daily max water level)
    surv_val = np.empty(len(med_h_wl)) #initiate an empty array
    for ii in range(h_wl.shape[0]):
        fromheightcalc, Pvaluecalc = calcWOO(h_wl[ii,:],woo_inun) # get correlation of elevation and Probability
        surv_val[ii] = np.interp(med_h_wl[ii],fromheightcalc,Pvaluecalc)
    
    # check surv_val, lower than LLWL should be 0
    bed_is = model_dfm.get_var('bl')
    surv_is = copy.copy(surv_val)
    for surv in range(len(surv_val)):
        if bed_is[surv] >= LLWL:
            surv_val[surv] = surv_is[surv]
        else:
            surv_val[surv] = 0
        
    ## Filter seedlings based on the surv_val
    try:
        if seedling_finalpos.size > 0:
            seedling_finalpos_2 = seedling_prob(seedling_finalpos, xzyz_cell_number,
                                        ugrid_all, surv_val)
            seedling_finalpos_filt = elim_seeds_surv(seedling_finalpos_2, xzyz_cell_number, ugrid_all, surv_val)
            seedling__finalpos = copy.copy(seedling_finalpos_filt)
    except:
        print('No seedlings establishment')
            
    ### 4. Convert from data point to raster environment 
    # 4.1 Create raster from the surv-val
    surv_val_raster = np.column_stack((xz,yz,surv_val))
    concave_path = os.path.join(MFON_Exchange, 'coupling'+str(ntime+1))
    if not os.path.exists(concave_path):
        os.makedirs(concave_path)
    surv_val_name = 'CH_surv_val'
    csv2ClippedRaster(concave_path, surv_val_raster, surv_val_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_Project)
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
    createPointSHP(read_data, age_coupling, concave_path, EPSG_Project) 
    
    # 6.2. Tile the Master Trees and Save as XLS Input File
    shp_source = os.path.join(concave_path, 'coupling'+str(ntime+1)+'.shp') # location of the source tree and the master tree shp
    folder_loc = dir_out # location of the master tiles
    file_tile = os.path.join(dir_out,'tile_*.shp' )
    save_tiled_trees = os.path.join(MFON_Trees,'coupling'+str(ntime+1)) # location to save the tiled shp
    if not os.path.exists(save_tiled_trees):
        os.makedirs(save_tiled_trees)
      
    # clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name, a0, b0, a137, b137)
    New_clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name)
        
    ### 8. Create the env raster and tile raster files
    save_tiled_env = os.path.join(MFON_Env, 'coupling'+str(ntime+1))
    if not os.path.exists(save_tiled_env):
        os.makedirs(save_tiled_env)
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    # calc_surv = '"0*(A==-999)"'
    # calc_sal = '"0*(A>-999)+35"'
    # val_no_data_surv = 0
    # val_no_data_sal = 60
    
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
        all_df.append(df)

    Concat_table = pd.concat(all_df)
    
    # 12.1. drop tick 0 year, only take 0.25
    Concat_table = Concat_table[Concat_table.tick > 0]
    run_is = 'Coupling_'+str(ntime+1) # change this with the real name
    
    # 12.2. use spatial in scipy to match the x,y of the mangroves and the age information.
    master_trees = gpd.read_file(os.path.join(concave_path+str('\\')+Path(concave_path).stem+'.shp'))
    age_coupling = calcAgeCoupling0(Concat_table, master_trees) # prepare the age_coupling for the Master Trees of Each Coupling procedure
    Concat_table['Age'] = age_coupling #update age after MesoFON run
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
    
### End Loop
#Finalize the running
model_dimr.finalize()  



