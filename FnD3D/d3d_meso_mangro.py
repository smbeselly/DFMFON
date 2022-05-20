# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 19:38:02 2021

@author: sbe002
"""

#%% Import Libraries
import numpy as np
import numpy.ma as ma
import pandas as pd
from sklearn.utils import resample
from scipy import spatial
import fiona
from pathlib import Path
from scipy.interpolate import interp1d
import glob
import geopandas as gpd
from xlrd import open_workbook
from xlutils.copy import copy
import os
from osgeo import gdal, gdalconst
import shutil
import fileinput
import re
import itertools
import sys
sys.path.append('D:/Git/d3d_meso/FnD3D')
from d3d_prep_raster import d3dPolySHP

import pandas as pd
from xlrd import open_workbook
import rasterio as rio
from affine import Affine

#%% Function to rearrange the cell number based on the xzw and yzw position

def create_xyzwCellNumber(xzw,yzw,mesh_face_x,mesh_face_y):
    """Create the matrix that contains the cell number based on the internal
    cell centre
    
    Parameters
    ---------
    xzw = x coord of centre cell of the internal elements/ flowcell
    yzw = y coord of centre cell of the internal elements/ flowcell
    
    xzw and yzw can be derived with this command
    
    xzw = model_dfm.get_var('xzw')    
    yzw = model_dfm.get_var('yzw')
    
    mesh_face_x = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_x')
    mesh_face_y = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_y')
    
    Returns
    ---------
    xyzw_cell_number matrix    
    shape: ndxi, 3    
    column: xzw, yzw, cell number
    """
    # xzw = model_dfm.get_var('xzw') #x coordinate of the center of gravity of the boxes
    # xzw = model_dfm.get_var('xzw') #y coordinate
    
    xyzw = np.block([[xzw],[yzw]]).T #matrix xzw,yzw
    mesh_face_xy = np.block([[mesh_face_x],[mesh_face_y]]).T #matrix mesh_face_x, mesh_face_y
    
    #create matrix: xzw,yzw,cell number
    xyzw_cell_number = np.zeros((len(xyzw),3)) #create zeros array
    xyzw_cell_number[:,0:2]=xyzw #broadcast xyzw array to :,0:2
    
    from scipy import spatial
    
    # assign the correct cell number to the matrix (cell number is python based)
    # if the cell number is read as 0, the real cell number in netCDF is 1
    #source: https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
    for row in range(len(xyzw)): #source: https://www.pluralsight.com/guides/numpy-arrays-iterating
        pt = [xyzw[row,:]]
        index = spatial.KDTree(mesh_face_xy).query(pt)
        xyzw_cell_number[row,2] = index[1]
    
    return xyzw_cell_number

def create_xzyzCellNumber(xz, yz, model_dfm, ugrid_all):
    xzyz = np.block([[xz],[yz]]).T 
    xzyz = xzyz[:model_dfm.get_var('ndxi')]
    
    df_xzyz = pd.DataFrame(xzyz)
    
    data_verts = ugrid_all.verts
    
    vegsp_index = np.zeros((len(xzyz),2))
    vegsp_index[:,0] = np.arange(0,len(xzyz))
    
    for row in range(len(data_verts[:,0,0])):
        aa = data_verts[row,:,:]
        
        ab_subset = df_xzyz[(df_xzyz[0] > min(aa[:,0])) & 
                                     (df_xzyz[0] < max(aa[:,0]))]
        ab_subset = ab_subset[(ab_subset[1] > min(aa[:,1])) & 
                                            (ab_subset[1] < max(aa[:,1]))]
        vegsp_index[row,1] = ab_subset.index.to_numpy()
        # vegsp_index[row,1] = 3*row
        
    #sort the numpy array by xzyz index
    sorted_vegsp_index = vegsp_index[vegsp_index[:,1].argsort()]
    #create xzyz_cell_number
    xzyz_cell_number = np.zeros((len(vegsp_index),3))
    xzyz_cell_number[:,0:2] = xzyz
    xzyz_cell_number[:,2] = sorted_vegsp_index[:,0]
    
    return xzyz_cell_number

#%% Function to rearrange the mesh_face_nodes that previously based on cell number
# now rearrange this based on the ndxi index -> create new matrix

def create_xyzwNodes(mesh_face_nodes,xyzw_cell_number):
    """Rearrange the mesh_face_nodes with the new one that follows the index
    of the xzw and yzw.
    
    So that we can directly use the node coordinates here to filter the mangrove 
    trees and place the calculated value in the same location as values read
    by the BMI
    
    Parameters
    ---------
    mesh_face_nodes = nodes number as in the domain (netCDF)
    read with -> 
    mesh_face_nodes = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_nodes')    
    xyzw_cell_number = matrix created from create_xyzwCellNumber
    
    Returns
    ---------
    xyzw_nodes    
    shape: ndxi, 4 or more column dependent on the mesh_face_nodes shape    
    column:all nodes number/ position
    """
    # mesh_face_nodes = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_face_nodes')
    
    # xyzw_nodes contains the nodes information based on the bmi-retrieved information.
    xyzw_nodes = ma.empty((ma.shape(mesh_face_nodes))) # create the empty False masked array 
    ma.set_fill_value(xyzw_nodes, -999)
    xyzw_nodes[:,:] = ma.masked # set an all True masked array first
    
    for row in range(len(xyzw_nodes)):
        xyzw_nodes[row] = mesh_face_nodes[xyzw_cell_number[row,2].astype(int)]
    
    return xyzw_nodes
    # now we have nodes information based on the xzw and yzw position in xyzw_nodes
    # if we want to check the correlating cell number of the xzw and yzw position,
    # one can check the xyzw_cell_number in column 2

# assume mangrove 

#%% Function to read trees simulation and do the drag coefficient calculation

def calcDragCoeff(x_range, y_range, cell_area, water_depth, trees_data):
    """Calculate the spatially varying drag coefficient.
    
    It is located in the cell centre
    
    
    This function is placed under the for-loop to calculate the drag coefficient
    for each cell.
    
    Currently, it only calculates Avicennia Marina. However an advancement
    can be made to accommodate another species. Target to upgrade the function
    to calculate Bruguiera (gymnorhiza) with its root knees, and Kandelia (candel)
    with the plank roots.
    
    Parameters
    ---------
    x_range = list of x min x max of the cell
    
    y_range = list of y min y range of the cell
    
    cell_area = area of the cell (m^2)
    
    water_depth = water depth (m)
    
    trees_data = reading of the tree simulation (read the GeoRefPosX, GeoRefPosY,
    dbh_cm)
    
    Returns
    ---------
    Cd_calc_weighted
    
    Calculate the weighted bulk drag coefficient of the cell, located in
    cell centre
    """
    ## If we want to add different species, we can use this:
# =============================================================================
#         import pandas as pd
#         data_url = 'http://bit.ly/2cLzoxH'
#         # read data from url as pandas dataframe
#         gapminder = pd.read_csv(data_url)
# 
#         searchstr = r'In(?!$)'
# 
#         data = gapminder[gapminder['country'].str.contains(searchstr)]
#   source : https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
# =============================================================================
    ## untuk test saja dummy value
    x_range = x_range
    y_range = y_range
    read_data = trees_data

    # subsetting pandas 
    #source: https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
    # and https://www.geeksforgeeks.org/selecting-rows-in-pandas-dataframe-based-on-conditions/
    read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                 (read_data['GeoRefPosX'] <= x_range[1])]
    read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                        (read_data_subset['GeoRefPosY'] <= y_range[1])]

    #retrieve cell area and water depth
    cell_area = cell_area # cell area based on BMI 
    water_depth = water_depth # water depth based on BMI

    # Calculate the CPRS and Number of Pneumatophore
    # Define the boundary (CPRS) as in Vovides, et al.,2016
    a_cprs = 1.45 
    b_cprs = 8.97
    cprs_avicennia = (a_cprs*read_data_subset['dbh_cm'])/(b_cprs+read_data_subset['dbh_cm']) #results in meter

    # N Pneumatophore as in van Maanen 
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen

    # Volume of the aerial roots
    N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-read_data_subset['dbh_cm']))))
    
    
    #if function to adjust the h_pneu for Avicennia
    # this equation is from Du, Qin, et al. 2021
    if water_depth < h_pneu:
        Vm_pneu = np.pi*d_pneu*water_depth/12
    else:
        Vm_pneu = np.pi*d_pneu*h_pneu/12

    Vm_pneu_total = Vm_pneu*N_pneu #m^3
    Am_pneu_total = d_pneu*h_pneu*N_pneu #m^2

    # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
    # height is adjusted from the water depth
    Vm_trunk = np.pi/4*(read_data_subset['dbh_cm']/100)**2*water_depth #m^3 
    Am_trunk = (read_data_subset['dbh_cm']/100)*water_depth #m^3 

    # sum of the obstacle volume
    Vm = Vm_pneu_total + Vm_trunk #m^3
    # sum of the obstacle area
    Am = Am_pneu_total + Am_trunk#m^2

    # Volume of the control water

    V = cell_area*water_depth

    # Characteristic Length (m)
    L = (V-Vm)/Am
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
    e_drag = 5 # set to 5m to obtain realistic value for Cd
    if L.size == 0:
        Cd_calc_weighted = Cd_no
    else:
        Cd_calc = Cd_no + (e_drag/L)
        Cd_calc_weighted = cprs_avicennia/ cprs_avicennia.sum() * Cd_calc
        Cd_calc_weighted = Cd_calc_weighted.sum()
    
    return Cd_calc_weighted

def calcWOO(waterlevel, failure):
    """Calculate probability of Windows of Opportunity.
    
    It calculates WoO based on the daily highest water level
    
    
    This function is placed under the for-loop to calculate the WoO Probability
    on each water level (s1) at cell center.
    
    It is based on Balke et.al.,2014
    Critical transitions in disturbance-driven ecosytems: identifying 
    Windows of Opportunity for recovery
          
    Parameters
    ---------
    waterlevel = daily highest water level (s1)
    
    failure = period free from inundation which needs to be exceeded for seedlings to anchor after stranding
       
    Returns
    ---------
    Pvalue
    
    WoO probability value of the waterlevel
    """
    
    h = waterlevel

    failure = failure # failure= period free from inundation which needs to be exceeded for seedlings to anchor after stranding
    diff = 0
    maxdur = np.zeros(len(h)) ## array for output: inundation free period
    maxel = np.zeros(len(h)) # array for output: water level which inundation free period is calculated for

    for i in range(len(h)):
        maxdur[i] = 0
        maxel[i] = h[i]
        if i < len(h)-1: # to prevent boundary problem when assess current to next
            diff = h[i] > h[i+1]
            if diff == True:
                q = i + 1
                while (h[i] > h[q]) and (q < len(h)-1):
                    maxdur[i] = maxdur[i] + 1
                    q += 1

                
    WoONa = np.block([[maxel],[maxdur]]).T
    WoO = WoONa[~np.isnan(WoONa).any(axis=1)]
    disperse = WoO[(WoO[:,1] >= 1)] #selects all days where diaspore strand all WoO >= 1 day no inundation,

    # convert disperse to dataframe to get density kernel
    disperse_df = pd.DataFrame(disperse)
    disperse_df.columns =['maxel','maxdur']
    
    # Calculate P Success
    establish = np.zeros(len(disperse))

    for i in range(len(disperse)):
        dispersal = disperse[i,1] > failure
        seedling = dispersal == True
        if seedling == True:
            establish[i] = 1
        else:
            establish[i] = 0

    WoOfin = np.column_stack((establish,disperse))

    ### probability of survival for different elevations, P success ####
    d1 = WoOfin[(WoOfin[:,0] == 1)]
    step = np.arange(np.ceil(np.sqrt(len(d1)))+1)
    ## sqrt of N,round up, define number of steps for probability analysis
    ## along height gradient
    
    Pvalue = np.zeros(int(len(step)))
    fromheight = np.zeros(int(len(step)))
    toheight = np.zeros(int(len(step)))
    PvalueNo = np.zeros(int(len(step)))
    startel = min(WoOfin[:,1])
    stepsize = (max(WoOfin[:,1]) - min(WoOfin[:,1]))/len(step)

    ## calculate size of steps probability calculation
    # for i in range(int(max(step))):
    for i in range(int(len(step))):
        Pselect2 = 0
        Pselect1 = 0
        Pselect2 = WoOfin[(WoOfin[:,1] >= (startel + (stepsize*i)))]
        Pselect1 = Pselect2[(Pselect2[:,1] <= (startel + (stepsize*(1+i))))]
        if Pselect1.size > 0:
            Pvalue[i] = sum(Pselect1[:,0])/len(Pselect1[:,0])
            PvalueNo[i] = sum(Pselect1[:,0])/len(Pselect1[:,0]) * sum(Pselect1[:,0])
        else:
            Pvalue[i] = 0
            PvalueNo[i] = 0
        # Pvalue[i] = sum(Pselect1[:,0])/len(Pselect1[:,0])
        # PvalueNo[i] = sum(Pselect1[:,0])/len(Pselect1[:,0]) * sum(Pselect1[:,0])
        fromheight[i] = startel + (stepsize * (i-1))
        toheight[i] = startel + (stepsize * (i))
    
    return fromheight, Pvalue

def calcAgeCoupling0(read_data, master_trees):
    """Calculate total age after initialization.
    
          
    Parameters
    ---------
    read_data = mangrove compilation data after MesoFON run
    
    master_trees = mangrove data from master trees
       
    Returns
    ---------
    age_coupling0
    
    total age after initialization
    
    Index is based on the read_data
    """
    read_data_np = read_data[['GeoRefPosX','GeoRefPosY']].to_numpy()
    master_trees_np = master_trees[['coord_x','coord_y']].to_numpy()
    # assign the correct cell number to the matrix (cell number is python based)
    # if the cell number is read as 0, the real cell number in netCDF is 1
    #source: https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
    age_read_data = np.empty((len(read_data),0))
    for row in range(len(read_data)): #source: https://www.pluralsight.com/guides/numpy-arrays-iterating
        pt = [read_data_np[row,:]]
        index = spatial.KDTree(master_trees_np).query(pt)
        age_read_data = np.append(age_read_data,master_trees['age'][int(index[1])])
        # age_read_data.append(master_trees['age'][int(index[1])])

    age_coupling0 = read_data['tick']+age_read_data
    age_coupling0 = age_coupling0.reset_index(drop=True)
    
    return age_coupling0

def createPointSHP(read_data, age_coupling, path_shp, EPSG_Project):
    """Create a Point Shapefile
    
          
    Parameters
    ---------
    read_data = mangrove compilation data after MesoFON run
    
    age_coupling = age data from age_coupling
    
    path_shp = path to save the shapefile
    
    EPSG_Project = str EPSG
       
    Returns
    ---------
    Point SHP
    
    """
    concave_path = path_shp
    data_point = [read_data['GeoRefPosX'], read_data['GeoRefPosY'], 
                  read_data['Height_cm']/100,age_coupling,read_data['dbh_cm']/200]
    data_point_header = ['coord_x','coord_y','height_m','age','rbh_m']
    data_point_df = pd.concat(data_point, axis=1, keys=data_point_header)
    # define schema
    schema = {
        'geometry':'Point',
        'properties':{'coord_x':'float',
                      'coord_y':'float',
                      'height_m':'float',
                      'age':'float',
                      'rbh_m':'float'}
    }

    #open a fiona object
    pointShp = fiona.open(concave_path+str('\\')+Path(concave_path).stem+'.shp', mode='w', driver='ESRI Shapefile',
              schema = schema, crs = 'EPSG:'+str(EPSG_Project))
    #iterate over each row in the dataframe and save record
    for index, row in data_point_df.iterrows():
        rowDict = {
            'geometry' : {'type':'Point',
                         'coordinates': (row.coord_x,row.coord_y)},
            'properties': {'coord_x' : row.coord_x,
                           'coord_y' : row.coord_y,
                           'height_m' : row.height_m,
                           'age' : row.age,
                           'rbh_m' : row.rbh_m},
        }
        pointShp.write(rowDict)
    #close fiona object
    pointShp.close()

  
def createXLSfromSHP(file_tile_trees, a0, b0, a137, b137, save_tiled_trees, species_name):
    """Create xls for data input to MesoFON
    
          
    Parameters
    ---------
    file_tile_trees = the tiled pt shp of mangrove trees
    
    a0, b0, a137, b137 = growth parameter, currently only one species
    
    save_tiled_trees = path to save the xls in the same location as the tiled trees
    
    species_name = name of the species, defined in the parameters cell
       
    Returns
    ---------
    XLS files
    
    """
    
    # calculate the h137
    d137 = np.arange(0,125.5,0.5)
    d0 = d137
    h137 = a137*d137**2+b137*d137+137
    h0 = a0*d0**2+b0*d0
    # calculate the correlation with interp1d
    f_rbh = interp1d(h137,d137, fill_value="extrapolate")
    f_0 = interp1d(h0,d0)

    for filepath in glob.iglob(file_tile_trees): # looping for all with trees affix
        # print(filepath)
        tile_0_read = gpd.read_file(filepath)
        posX = tile_0_read['coord_x']
        posY = tile_0_read['coord_y']
        height_m = tile_0_read['height'] #change for independent files or standardized the shp
        id_id = np.arange(len(posX))
        speciesName = np.ones(len(posX))*1 #if only one species is recorded, however it is better to place this in shapefile
        types_species = id_id+1 # this variable starts from 1 to N+1
        # shiftedBelowPos = np.ones(len(posX))*1
        age = tile_0_read['age']
        # age = np.ones(len(posX))*1 # should be derived from shapefile #TO DO will be added later
        # rbh_m = (0.0015*height_m**2)+(0.0015*height_m) #dummy it uses equation in test_trial
        # use the new dbh-height relationship  
        if height_m.size != 0:
            height_m = height_m.values  # in metre
            rbh_m = np.where(height_m < 1.37, f_0(height_m*100)/100, f_rbh(height_m*100)/100)
        else:
            rbh_m = tile_0_read['height']

        # Create Panda Dataframe with Header
        df = pd.DataFrame({'id': id_id,
                           'speciesName' : speciesName,
                           'type' : types_species,
                           'posX' : posX,
                           'posY' : posY,
                           'shiftedPosX' : posX,
                           'shiftedPosY' : posY,
                           'shiftedBelowPosX' : speciesName,
                           'shiftedBelowPosY' : speciesName,
                           'age' : age,
                           'rbh_m' : rbh_m,
                           'newRbh_m' : rbh_m,
                           'height_m' : height_m,
                           'newHeight_m' : height_m})
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        xls_name = Path(filepath).stem+'_input'+'.xls'
        xls_loc = save_tiled_trees + xls_name
        
        # Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
        # the default header and skip one row to allow us to insert a user defined
        # header.
        startrowis = 4
        
        with pd.ExcelWriter(xls_loc) as writer:
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=startrowis, header=True)
        
        rb = open_workbook(xls_loc)
        wb = copy(rb)
        # Insert the value in worksheet
        s = wb.get_sheet(0)
        s.write(0,0,'def')
        s.write(1,0, 1)
        s.write(1,1, species_name)
        s.write(2,0, '/def')
        s.write(3,0, 'values')
            # position_last = start the data + 1 (since it contains header) 
        #  + length of data
        position_last = startrowis + 1 + len(posX)
        s.write(position_last,0, '/values')
        wb.save(xls_loc)

       
def createRaster4MesoFON(concave_path, gdal_calc_path, no_data_val, save_tiled_env, calc_sal, val_no_data_sal):
    """Create tiled surv and sal rasters for each coupling 
    
    It will create folder coupling+str(ntime+1) and copy the rasters to the folder
          
    Parameters
    ---------
    concave_path = path to the master tile (calculated surv)
    
    gdal_calc_path = path for gdal (it is manually referred)
    
    no_data_val = taken from the params cell -999.0
    
    save_tiled_env = path to save the tiled_env
    
    calc_sal = algorithm to create idealized salinity
    
    val_no_data_sal = value of no data for salinity
    
       
    Returns
    ---------
    tiled raster files
    
    """
    dir_out = concave_path
   
    def changeNoData(datatif,value):
        maskfile = gdal.Open(datatif, gdalconst.GA_Update)
        maskraster = maskfile.ReadAsArray()
        maskraster = np.where((maskraster > -999), maskraster, value ) # yg lebih dari -999 diganti jadi 0
        maskband = maskfile.GetRasterBand(1)
        maskband.WriteArray( maskraster )
        maskband.FlushCache()
        
    for filepatd in glob.iglob(os.path.join(dir_out,'tile_*.tif')):
        # print(filepatd)
        ori_raster = filepatd
        tif_name = Path(filepatd).stem
        # surv_raster = os.path.join(save_tiled_env, tif_name+'_surv_.tif')
        surv_raster = filepatd # surv file has been calculated
        sal_raster = os.path.join(save_tiled_env, tif_name+'_sal_.tif')
        # define the syntax
        # command_surv = 'python {gdal_calc_path} -A {ori_raster} --outfile={surv_raster} --calc={calc_surv} --NoDataValue={no_data_val}'
        command_sal = 'python {gdal_calc_path} -A {ori_raster} --outfile={sal_raster} --calc={calc_sal} --NoDataValue={no_data_val}'
        # calculate
        # os.system(command_surv.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
        #                               surv_raster=surv_raster, calc_surv=calc_surv, no_data_val=no_data_val))
        os.system(command_sal.format(gdal_calc_path=gdal_calc_path, ori_raster=ori_raster, \
                                     sal_raster=sal_raster, calc_sal=calc_sal, no_data_val=no_data_val))
        # change the nodata value into 0 for surv and 60 for sal
        # changeNoData(surv_raster,val_no_data_surv) #since the surv raster is calculated from WoO 
        changeNoData(sal_raster,val_no_data_sal)
        # copy and rename the new raster for surv
        raster_surv = surv_raster
        raster_surv_name = Path(raster_surv).stem
        target_ras_surv0 = os.path.join(save_tiled_env,raster_surv_name+'0.tif')
        target_ras_surv1 = os.path.join(save_tiled_env,raster_surv_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_surv, target_ras_surv0)
        shutil.copyfile(raster_surv, target_ras_surv1)
        # copy and rename the new raster for sal
        raster_sal = sal_raster
        raster_sal_name = Path(raster_sal).stem
        target_ras_sal0 = os.path.join(save_tiled_env,raster_sal_name+'0.tif')
        target_ras_sal1 = os.path.join(save_tiled_env,raster_sal_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_sal, target_ras_sal0)
        shutil.copyfile(raster_sal, target_ras_sal1)
        
def modifyParamMesoFON(MFON_HOME, MFON_Exchange, save_tiled_env, save_tiled_trees, 
                       Surv_Source, Sal_Source, Excel_Source, ntime):
    """Modify unrolledParamFile.txt, batch_params.xml, and parameters.xml
    
    It modifies the path to the surv, sal, and excel file
          
    Parameters
    ---------
    MFON_HOME = path to MFON_HOME
    
    MFON_Exchange = path to MFON_Exchange
    
    save_tiled_env =  path to the tiled_env
    
    save_tiled_trees = path to the tiled_trees
    
    Surv_Source = str of surv location in payload.jar
      
    Sal_Source = str of sal location in payload.jar
    
    Excel_Source = str of excel location in payload.jar
    
    ntime = ntime

    Returns
    ---------
    updated parameters file (unrolledParamFile, batch_params, and parameters)
    
    """
    # source = https://www.delftstack.com/howto/python/python-replace-line-in-file/
    def replacement(file, previousw, nextw):
       for line in fileinput.input(file, inplace=1):
           line = line.replace(previousw, nextw)
           sys.stdout.write(line)
    # The Looping to replace the variables in the params files
    for filepatf in glob.iglob(os.path.join(MFON_Exchange,'Initialization','tile_*.tif')):
        # Init_Rasters = save_tiled_env
        # Init_Trees = save_tiled_trees
        # Create Update Path
        # Sal_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_sal_').replace("\\",'/')
        # Surv_Update = os.path.join(Init_Rasters,Path(filepatf).stem+'_surv_').replace("\\",'/')
        # Excel_Update = os.path.join(Init_Trees,Path(filepatf).stem+'_trees_input.xls').replace("\\",'/')
        Sal_Update = 'coupling'+str(ntime+1)
        Surv_Update = 'coupling'+str(ntime+1)
        Excel_Update = 'coupling'+str(ntime+1)
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
    
    return Surv_Update, Sal_Update, Excel_Update

def calcDragInLoop(xyzw_cell_number, model_dfm, xyzw_nodes, xk, yk, read_data):
    # drag_coeff = [] #should explore whether as list or as array list
    drag_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))

    for row in range(len(xyzw_cell_number)):
        ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
        hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
        # find the position based on the cell number
        position = xyzw_cell_number[row,2].astype(int)
        
        nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
        nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
        # Find the min max of each x,y coordinate
        # create the list of x_min-x_max and y_min-y_max
        x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
        y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
        
        # subsetting pandas 
        read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                     (read_data['GeoRefPosX'] <= x_range[1])]
        read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                            (read_data_subset['GeoRefPosY'] <= y_range[1])]
        # TODO check cell_area and water_depth
        cell_area = ba[row] # check this later 
        water_depth = hs[row]# adjust
        trees_data = read_data_subset
        # calculate the drag coefficient (currently only for Avicennia marina)
        cd_veg = calcDragCoeff(x_range, y_range, cell_area, water_depth, trees_data)
        # append the calculated cd_veg to the drag_coeff list
        # drag_coeff.append(cd_veg)
        drag_coeff = np.append(drag_coeff, cd_veg)
    
    return drag_coeff

def csv2ClippedRaster(concave_path, surv_val_raster, concave_name, x_res, y_res, no_data_val, affix, dir_out, EPSG_coord):
    """New csv to raster and clipped raster
    
    As the old os.system approach in d3dCSV2ClippedRaster doesn't work
          
    Parameters
    ---------
    concave_path = path to coupling+str(ntime+1))
    
    surv_val_raster = surv_val_raster
    
    x_res =  x_res
    
    y_res = y_res
    
    no_data_val = no_data_val
      
    affix = the end of the file name designed as '_clipped'
    
    dir_out = path to the CH_.shp in initialization
    
    EPSG_coord = EPSG_coord

    Returns
    ---------
    updated parameters file (unrolledParamFile, batch_params, and parameters)
    
    """
    os.chdir(concave_path)
    
    matrix = surv_val_raster
    
    np.savetxt(str(concave_path)+str('\\')+str(concave_name)+'.csv', matrix, fmt="%f", delimiter=",", header='Lon,Lat,Alt',comments='')
    
    # Create .vrt file based on the created .csv
    lines = ['<OGRVRTDataSource>', '<OGRVRTLayer name='+'"'+str(concave_name)+'"'+'>',
             '<SrcDataSource>'+str(concave_name)+".csv"+'</SrcDataSource>',
             '<GeometryType>wkbPoint</GeometryType>',
             '<LayerSRS>EPSG:'+str(EPSG_coord)+'</LayerSRS>',
             '<GeometryField separator=" " encoding="PointFromColumns" x="Lon" y="Lat" z="Alt"/>',
             '</OGRVRTLayer>',
             '</OGRVRTDataSource>']
    with open(str(concave_path)+str('\\')+str(concave_name)+'.vrt', 'w') as f:
        f.write('\n'.join(lines))
     
    # x_min = np.min(surv_val_raster[:,0])
    # x_max = np.max(surv_val_raster[:,0])
    # y_min = np.min(surv_val_raster[:,1])
    # y_max = np.max(surv_val_raster[:,1])
    # x_width = (x_max-x_min)/x_res
    # y_height = (y_max-y_min)/y_res
    template = os.path.join(dir_out,'CH_bathy__clipped.tif')
    ds = gdal.Open(template)
    gt = ds.GetGeoTransform()
    res = gt[1]
    x_width = ds.RasterXSize
    y_height = ds.RasterYSize
    x_min =  gt[0]
    x_max =  x_min + res* x_width
    y_min = gt[3]
    y_max = y_min + gt[5] * y_height

    lin_raster = gdal.Grid(os.path.join(concave_path,concave_name+'.tif'), 
                            os.path.join(concave_path,concave_name+'.vrt'),
                            algorithm = "linear", noData = no_data_val ,zfield = 'Alt',
                    outputBounds = [x_min,y_max,x_max,y_min],
                    width = x_width, height = y_height)
    lin_raster = None

    ras_clip = str(concave_name)+affix+'.tif'
    cut_call = os.path.join(dir_out,'CH_.shp')

    warpp = gdal.Warp(os.path.join(concave_path,ras_clip), 
                      os.path.join(concave_path,concave_name+'.tif'), cutlineDSName = cut_call, 
                      cropToCutline = True, dstNodata=no_data_val)

    warpp = None
    
def d3dNewRaster2Tiles(ras_clip, out_path, tile_size_x, tile_size_y, dir_out):
    output_filename = "tile_"
    ds = gdal.Open(ras_clip)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    # stats = band.GetStatistics(True,True) # results = min, max, mean, StdDev
    # stats[1]-stats[0] = max-min
    xsize = band.XSize
    ysize = band.YSize
    # get coordinates of upper left corner
    # xmin = gt[0]
    # ymax = gt[3]
    res = gt[1]

    # close dataset
    ds = None
    # determine total length of raster
    # xlen = res * ds.RasterXSize # convert the size in meter to number of pixels
    # ylen = res * ds.RasterYSize

    # size of a single tile
    xsize_tile = round(tile_size_x/res) #num of pixels in tile_size_x
    ysize_tile = round(tile_size_y/res) ##num of pixels in tile_size_y

    # Tile the raster domain as the prefered tile size in meter
    for i in range(0, xsize, xsize_tile):
        for j in range(0, ysize, ysize_tile):
            print (i,j)
            prep_out = str(out_path) +str('\\')+ str(output_filename) + str(i) + "_" + str(j) + ".tif"  
            if os.path.isfile(os.path.join(dir_out,Path(prep_out).name)) == True:
                translate_tile = gdal.Translate(prep_out, ras_clip, srcWin = [i,j,xsize_tile,ysize_tile])
                translate_tile = None
            else:
                # below is to delete the grid with all nodata value
                # del(dss,bands,statss)
                prep_del = prep_out
                del(prep_out)
                os.remove(str(prep_del))
                
def Sald3dNewRaster2Tiles(ras_clip, out_path, tile_size_x, tile_size_y, dir_out, output_filename):
    ds = gdal.Open(ras_clip)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    # stats = band.GetStatistics(True,True) # results = min, max, mean, StdDev
    # stats[1]-stats[0] = max-min
    xsize = band.XSize
    ysize = band.YSize
    # get coordinates of upper left corner
    # xmin = gt[0]
    # ymax = gt[3]
    res = gt[1]

    # close dataset
    ds = None
    # determine total length of raster
    # xlen = res * ds.RasterXSize # convert the size in meter to number of pixels
    # ylen = res * ds.RasterYSize

    # size of a single tile
    xsize_tile = round(tile_size_x/res) #num of pixels in tile_size_x
    ysize_tile = round(tile_size_y/res) ##num of pixels in tile_size_y

    # Tile the raster domain as the prefered tile size in meter
    for i in range(0, xsize, xsize_tile):
        for j in range(0, ysize, ysize_tile):
            print (i,j)
            prep_out = str(out_path) +str('\\')+ str(output_filename) + str(i) + "_" + str(j) + ".tif" 
            check_dir = str(out_path) +str('\\')+ 'tile_' + str(i) + "_" + str(j) + ".tif" 
            if os.path.isfile(os.path.join(dir_out,Path(check_dir).name)) == True:
                translate_tile = gdal.Translate(prep_out, ras_clip, srcWin = [i,j,xsize_tile,ysize_tile])
                translate_tile = None
            else:
                # below is to delete the grid with all nodata value
                # del(dss,bands,statss)
                prep_del = prep_out
                del(prep_out)
                os.remove(str(prep_del))
                
def clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name, a0, b0, a137, b137):
    """Clip the Master Trees and Create XLS from the GPD read
    
          
    Parameters
    ---------
    file_tile = path to the tile for looping
    
    save_tiled_trees = path to save the tiled trees and xls
    
    shp_source =  master trees
    
    species_name = species name
    
    a0, b0, a137, b137 = Avicennia Marina Parameters

    Returns
    ---------
    XLS files for MesoFON run
    
    """

    for filepath in glob.iglob(file_tile):
        # print(filepath)
        # save_name = Path(filepath).stem+'_trees'+'.shp'
        # save_loc = os.path.join(folder_loc,save_name)
        # save_loc = os.path.join(save_tiled_trees,save_name)
    # =============================================================================
    #     command_ogr = 'ogr2ogr -clipsrc {filepath} {save_loc} {shp_source} -f "ESRI Shapefile"'
    #     os.system(command_ogr.format(filepath=filepath, save_loc=save_loc, 
    #                                   shp_source=shp_source))
    # =============================================================================
        
       
        d137 = np.arange(0,125.5,0.5)
        d0 = d137
        h137 = a137*d137**2+b137*d137+137
        h0 = a0*d0**2+b0*d0
        # calculate the correlation with interp1d
        f_dbh = interp1d(h137,d137, fill_value="extrapolate")
        f_0 = interp1d(h0,d0)
        # calculate the clipping with geopandas
        # gp_point = shp_source
        gp_point= gpd.read_file(shp_source)
        # your_clip = os.path.join(r"D:\Git\d3d_meso\Model-Exchange\Initialization\tile_0_20.shp")
        gp_clip= gpd.read_file(filepath)
    
        tree_point = gp_point.clip(gp_clip)
        # tree_point.to_file(save_loc)
        
        # create the XLS file based on the clipped shp
        posX = tree_point['coord_x']
        posY = tree_point['coord_y']
        height_m = tree_point['height_m'] #change for independent files or standardized the shp
        id_id = np.arange(len(posX))
        speciesName = np.ones(len(posX))*1 #if only one species is recorded, however it is better to place this in shapefile
        types_species = id_id+1 # this variable starts from 1 to N+1
        # shiftedBelowPos = np.ones(len(posX))*1
        age = tree_point['age']
        # age = np.ones(len(posX))*1 # should be derived from shapefile #TO DO will be added later
        # rbh_m = (0.0015*height_m**2)+(0.0015*height_m) #dummy it uses equation in test_trial
        # use the new dbh-height relationship  
        if height_m.size != 0:
            # height_m = height_m.values  # in metre
            rbh_m = np.where(height_m < 1.37, f_0(height_m*100)/100/2, f_dbh(height_m*100)/100/2)
            # (height_m*100)/100/2 karena mengonversi dari dbh ke rbh
            # rbh_m = np.where(height_m < 1.37, f_0(height_m*50)/50, f_rbh(height_m*50)/50)
        else:
            rbh_m = tree_point['height_m']
    
        # Create Panda Dataframe with Header
        df = pd.DataFrame({'id': id_id,
                           'speciesName' : speciesName,
                           'type' : types_species,
                           'posX' : posX,
                           'posY' : posY,
                           'shiftedPosX' : posX,
                           'shiftedPosY' : posY,
                           'shiftedBelowPosX' : speciesName,
                           'shiftedBelowPosY' : speciesName,
                           'age' : age,
                           'rbh_m' : rbh_m,
                           'newRbh_m' : rbh_m,
                           'height_m' : height_m,
                           'newHeight_m' : height_m})
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        xls_name = Path(filepath).stem+'_trees_input'+'.xls'
        xls_loc = os.path.join(save_tiled_trees , xls_name)
        
        # Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
        # the default header and skip one row to allow us to insert a user defined
        # header.
        startrowis = 4
        
        with pd.ExcelWriter(xls_loc) as writer:
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=startrowis, header=True)
        
        rb = open_workbook(xls_loc)
        wb = copy(rb)
        # Insert the value in worksheet
        s = wb.get_sheet(0)
        s.write(0,0,'def')
        s.write(1,0, 1)
        s.write(1,1, species_name)
        s.write(2,0, '/def')
        s.write(3,0, 'values')
            # position_last = start the data + 1 (since it contains header) 
        #  + length of data
        position_last = startrowis + 1 + len(posX)
        s.write(position_last,0, '/values')
        wb.save(xls_loc)

def _new_func_createRaster4MesoFON(concave_path,save_tiled_env, no_data_val, EPSG_Project, val_no_data_sal):
    dir_coupling = concave_path
    
    def changeNoData(datatif,value):
        maskfile = gdal.Open(datatif, gdalconst.GA_Update)
        maskraster = maskfile.ReadAsArray()
        maskraster = np.where((maskraster > -999), maskraster, value ) # yg lebih dari -999 diganti jadi 0
        maskband = maskfile.GetRasterBand(1)
        maskband.WriteArray( maskraster )
        maskband.FlushCache()
        
    for filepatd in glob.iglob(os.path.join(dir_coupling,'tile_*.tif')):
        ori_raster = filepatd
        tif_name = Path(filepatd).stem
        surv_raster = os.path.join(dir_coupling, tif_name+'.tif')
        # surv_raster = filepatd # surv file has been calculated
        sal_raster = os.path.join(save_tiled_env, tif_name+'_sal_.tif')
        # calc raster with rasterio and gdal
         
        ds = gdal.Open(filepatd, gdal.GA_ReadOnly)
        GT_input = ds.GetGeoTransform()
        afn = Affine.from_gdal(*GT_input)
        band = ds.GetRasterBand(1)
        raster = band.ReadAsArray()
        # define the ideal value
        raster[raster > no_data_val] = 35
        with rio.open(sal_raster, 'w', driver='GTiff', 
                           height=raster.shape[0], width=raster.shape[1],
                           count=1, dtype=np.float64, crs='EPSG:'+str(EPSG_Project),
                           transform=afn,) as dest_file: #nodata=no_data_val
            dest_file.write(raster, 1)
            dest_file.close()
        
        
        # change the nodata value into 0 for surv and 60 for sal
        # changeNoData(surv_raster,val_no_data_surv) #since the surv raster is calculated from WoO 
        changeNoData(sal_raster,val_no_data_sal)
        # copy and rename the new raster for surv
        raster_surv = surv_raster
        raster_surv_name = Path(raster_surv).stem+'_surv_'
        target_ras_surv0 = os.path.join(save_tiled_env,raster_surv_name+'0.tif')
        target_ras_surv1 = os.path.join(save_tiled_env,raster_surv_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_surv, os.path.join(save_tiled_env,raster_surv_name+'.tif'))
        shutil.copyfile(raster_surv, target_ras_surv0)
        shutil.copyfile(raster_surv, target_ras_surv1)
        # copy and rename the new raster for sal
        raster_sal = sal_raster
        raster_sal_name = Path(raster_sal).stem
        target_ras_sal0 = os.path.join(save_tiled_env,raster_sal_name+'0.tif')
        target_ras_sal1 = os.path.join(save_tiled_env,raster_sal_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_sal, target_ras_sal0)
        shutil.copyfile(raster_sal, target_ras_sal1)

### buat list untuk read data subset
def list_subset(xyzw_cell_number, index_veg_cel, xyzw_nodes, xk, yk, read_data):
    # this function is to predefine the subset before calculate drag
        
    list_read_subset = []
    for row in range(len(xyzw_cell_number)):
        if index_veg_cel[row] == 0:
            read_data_subset = np.nan
        else:
            position = xyzw_cell_number[row,2].astype(int)
            
            nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
            nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
            # Find the min max of each x,y coordinate
            # create the list of x_min-x_max and y_min-y_max
            x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
            y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
            
            # subsetting pandas 
            read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                         (read_data['GeoRefPosX'] <= x_range[1])]
            read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                                (read_data_subset['GeoRefPosY'] <= y_range[1])]
        
        list_read_subset.append(read_data_subset)
        
    return list_read_subset

def list_subsetCdveg(ugrid_all, xzyz_cell_number, index_veg_cel, read_data):
    # this function is to predefine the subset before calculate drag
    data_verts = ugrid_all.verts    
    list_read_subset = []
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 0:
            read_data_subset = np.nan
        else:
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                         (read_data['GeoRefPosX'] <= max(aa[:,0]))]
            read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                                (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
        
        list_read_subset.append(read_data_subset)
        
    return list_read_subset

### drag predictor formula
def initCalcDraginLoop(xyzw_cell_number, xyzw_nodes, xk, yk, read_data, model_dfm, index_veg_cel):
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
    # Parameters of the CPRS and Number of Pneumatophore as in Vovides, et al.,2016
    # a_cprs = 1.45 
    # b_cprs = 8.97
    # N Pneumatophore as in van Maanen 
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen
    ## Create an array of vegetated cell
    # index_veg_cel = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    # for row in range(len(xyzw_cell_number)):
    #     position = xyzw_cell_number[row,2].astype(int)
        
    #     nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
    #     nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
    #     # Find the min max of each x,y coordinate
    #     # create the list of x_min-x_max and y_min-y_max
    #     x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
    #     y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
        
    #     # subsetting pandas 
    #     read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
    #                                  (read_data['GeoRefPosX'] <= x_range[1])]
    #     read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
    #                                         (read_data_subset['GeoRefPosY'] <= y_range[1])]
    #     if read_data_subset.shape[0] == 0: #if 0 then skip
    #         index_is = 0
    #     else:
    #         index_is = 1
    #     index_veg_cel = np.append(index_veg_cel, index_is)
    # the boundary flow nodes are located at the end of array
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    
    drag_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    for row in range(len(xyzw_cell_number)):
        if index_veg_cel[row] == 0:
            Cd_calc_bare = 0.005
            drag_coeff = np.append(drag_coeff, Cd_calc_bare)
        else:
            # calculate cell_area and water_depth
            cell_area = ba[row] # sudah sesuai karena boundary ada di array paling akhir
            water_depth = hs[row]# sama seperti di atas
            # sometimes there is a condition where the cell is dry
            # therefore, if it is dry I set the drag as Cd_no
            if water_depth == 0:
                Cd_calc_bare = 0.005
                drag_coeff = np.append(drag_coeff, Cd_calc_bare)
            else:
                # print(row)
                # find the position based on the cell number
                position = xyzw_cell_number[row,2].astype(int)
                
                nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
                nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
                # Find the min max of each x,y coordinate
                # create the list of x_min-x_max and y_min-y_max
                x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
                y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
                
                # subsetting pandas 
                read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                             (read_data['GeoRefPosX'] <= x_range[1])]
                read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                                    (read_data_subset['GeoRefPosY'] <= y_range[1])]
               
                # calculate the drag coefficient (currently only for Avicennia marina)
                # cd_veg = calcDragCoeff(x_range, y_range, cell_area, water_depth, trees_data)
            
                # Calculate the CPRS and Number of Pneumatophore
                # Define the boundary (CPRS) as in Vovides, et al.,2016          
                # cprs_avicennia = (a_cprs*read_data_subset['dbh_cm'])/(b_cprs+read_data_subset['dbh_cm']) #results in meter
            
                # Volume of the aerial roots
                N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-read_data_subset['dbh_cm']))))
                          
                # if function to adjust the h_pneu for Avicennia
                # this equation is from Du, Qin, et al. 2021
                if water_depth < h_pneu:
                    Vm_pneu = np.pi*d_pneu*water_depth/12
                else:
                    Vm_pneu = np.pi*d_pneu*h_pneu/12
            
                Vm_pneu_total = Vm_pneu*N_pneu #m^3
                Am_pneu_total = d_pneu*h_pneu*N_pneu #m^2
            
                # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
                # height is adjusted from the water depth
                Vm_trunk = np.pi/4*(read_data_subset['dbh_cm']/100)**2*water_depth #m^3 
                Am_trunk = (read_data_subset['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am = Am_pneu_total + Am_trunk#m^2
            
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V-Vm)/Am
                # Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
                e_drag = 5 # set to 5m to obtain realistic value for Cd
            
                Cd_calc = (Cd_no + (e_drag/L)).sum()
                # TODO this is just for reminder of the commented script
                # is it really needed to calculate the weighted drag?
                # Cd_calc_weighted = Cd_calc*(Vm.sum()/V) + Cd_no * ((V-Vm.sum())/V)
                
                # append the calculated cd_veg to the drag_coeff list
                # drag_coeff.append(cd_veg)
                drag_coeff = np.append(drag_coeff, Cd_calc)
        
    return drag_coeff

def initCalcDraginLoopCdveg(xzyz_cell_number, model_dfm, ugrid_all, index_veg_cel, read_data):
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation 
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    
    data_verts = ugrid_all.verts
    
    drag_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 0:
            Cd_calc_bare = 0.005
            drag_coeff = np.append(drag_coeff, Cd_calc_bare)
        else:
            # calculate cell_area and water_depth
            cell_area = ba[row] # sudah sesuai karena boundary ada di array paling akhir
            water_depth = hs[row]# sama seperti di atas
            # sometimes there is a condition where the cell is dry
            # therefore, if it is dry I set the drag as Cd_no
            if water_depth == 0:
                Cd_calc_bare = 0.005
                drag_coeff = np.append(drag_coeff, Cd_calc_bare)
            else:
                # print(row)
                # find the position based on the cell number
                position = xzyz_cell_number[row,2].astype(int)
                aa = data_verts[position,:,:]
                # subsetting pandas 
                read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                             (read_data['GeoRefPosX'] <= max(aa[:,0]))]
                read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                                    (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
               
                # Volume of the aerial roots
                N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-read_data_subset['dbh_cm']))))
                          
                # if function to adjust the h_pneu for Avicennia
                # this equation is from Du, Qin, et al. 2021
                if water_depth < h_pneu:
                    Vm_pneu = np.pi*d_pneu*water_depth/12
                else:
                    Vm_pneu = np.pi*d_pneu*h_pneu/12
            
                Vm_pneu_total = Vm_pneu*N_pneu #m^3
                Am_pneu_total = d_pneu*h_pneu*N_pneu #m^2
            
                # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
                # height is adjusted from the water depth
                Vm_trunk = np.pi/4*(read_data_subset['dbh_cm']/100)**2*water_depth #m^3 
                Am_trunk = (read_data_subset['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am = Am_pneu_total + Am_trunk#m^2
            
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V-Vm)/Am
                # Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
                e_drag = 5 # set to 5m to obtain realistic value for Cd
            
                Cd_calc = (Cd_no + (e_drag/L)).sum()
                # TODO this is just for reminder of the commented script
                # is it really needed to calculate the weighted drag?
                # Cd_calc_weighted = Cd_calc*(Vm.sum()/V) + Cd_no * ((V-Vm.sum())/V)
                
                # append the calculated cd_veg to the drag_coeff list
                # drag_coeff.append(cd_veg)
                drag_coeff = np.append(drag_coeff, Cd_calc)
        
    return drag_coeff

def newCalcDraginLoop(xyzw_cell_number, read_data, model_dfm, index_veg_cel, list_read_subset):
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
    # Parameters of the CPRS and Number of Pneumatophore as in Vovides, et al.,2016
    # a_cprs = 1.45 
    # b_cprs = 8.97
    # N Pneumatophore as in van Maanen 
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen
    ## Create an array of vegetated cell
    # index_veg_cel = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    # for row in range(len(xyzw_cell_number)):
    #     position = xyzw_cell_number[row,2].astype(int)
        
    #     nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
    #     nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
    #     # Find the min max of each x,y coordinate
    #     # create the list of x_min-x_max and y_min-y_max
    #     x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
    #     y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
        
    #     # subsetting pandas 
    #     read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
    #                                  (read_data['GeoRefPosX'] <= x_range[1])]
    #     read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
    #                                         (read_data_subset['GeoRefPosY'] <= y_range[1])]
    #     if read_data_subset.shape[0] == 0: #if 0 then skip
    #         index_is = 0
    #     else:
    #         index_is = 1
    #     index_veg_cel = np.append(index_veg_cel, index_is)
    # the boundary flow nodes are located at the end of array
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    
    drag_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    rnveg_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    diaveg_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    stemheight_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    
    for row in range(len(xyzw_cell_number)):
        if index_veg_cel[row] == 0:
            Cd_calc_bare = 0.005
            drag_coeff = np.append(drag_coeff, Cd_calc_bare)
            rnveg_coeff = np.append(rnveg_coeff, 0)
            diaveg_coeff = np.append(diaveg_coeff, 0)
            stemheight_coeff = np.append(stemheight_coeff, 0)
        else:
            # calculate cell_area and water_depth
            cell_area = ba[row] # sudah sesuai karena boundary ada di array paling akhir
            water_depth = hs[row]# sama seperti di atas
            # sometimes there is a condition where the cell is dry
            # therefore, if it is dry I set the drag as Cd_no
            if water_depth == 0:
                Cd_calc_bare = 0.005
                drag_coeff = np.append(drag_coeff, Cd_calc_bare)
                rnveg_coeff = np.append(rnveg_coeff, 0)
                diaveg_coeff = np.append(diaveg_coeff, 0)
                stemheight_coeff = np.append(stemheight_coeff, 0)
            else:
                # print(row)
               # find the position based on the cell number
# =============================================================================
#                 position = xyzw_cell_number[row,2].astype(int)
#                 
#                 nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
#                 nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
#                 # Find the min max of each x,y coordinate
#                 # create the list of x_min-x_max and y_min-y_max
#                 x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
#                 y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
#                 
#                 # subsetting pandas 
#                 read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
#                                              (read_data['GeoRefPosX'] <= x_range[1])]
#                 read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
#                                                     (read_data_subset['GeoRefPosY'] <= y_range[1])]
# =============================================================================
                
                read_data_subset = list_read_subset[row]
                
    
                # calculate the drag coefficient (currently only for Avicennia marina)
                # cd_veg = calcDragCoeff(x_range, y_range, cell_area, water_depth, trees_data)
            
                # Calculate the CPRS and Number of Pneumatophore
                # Define the boundary (CPRS) as in Vovides, et al.,2016          
                # cprs_avicennia = (a_cprs*read_data_subset['dbh_cm'])/(b_cprs+read_data_subset['dbh_cm']) #results in meter
            
                # Volume of the aerial roots
                N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-read_data_subset['dbh_cm']))))
                          
                # if function to adjust the h_pneu for Avicennia
                # this equation is from Du, Qin, et al. 2021
                if water_depth < h_pneu:
                    Vm_pneu = np.pi*d_pneu*water_depth/12
                else:
                    Vm_pneu = np.pi*d_pneu*h_pneu/12
            
                Vm_pneu_total = Vm_pneu*N_pneu #m^3
                Am_pneu_total = d_pneu*h_pneu*N_pneu #m^2
            
                # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
                # height is adjusted from the water depth
                Vm_trunk = np.pi/4*(read_data_subset['dbh_cm']/100)**2*water_depth #m^3 
                Am_trunk = (read_data_subset['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am = Am_pneu_total + Am_trunk#m^2
            
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V-Vm)/Am
                # Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
                e_drag = 5 # set to 5m to obtain realistic value for Cd
            
                Cd_calc = (Cd_no + (e_drag/L)).sum()
                # TODO this is just for reminder of the commented script
                # is it really needed to calculate the weighted drag?
                # Cd_calc_weighted = Cd_calc*(Vm.sum()/V) + Cd_no * ((V-Vm.sum())/V)
                
                # append the calculated cd_veg to the drag_coeff list
                # drag_coeff.append(cd_veg)
                drag_coeff = np.append(drag_coeff, Cd_calc)
                
                # calculate rnveg, diaveg, and stemheight
                rnveg_calc = read_data_subset.Height_cm.notnull().sum()
                diaveg_calc = read_data_subset.dbh_cm.quantile(0.7)/100
                stemheight_calc = read_data_subset.Height_cm.quantile(0.7)/100
                
                rnveg_coeff = np.append(rnveg_coeff, rnveg_calc)
                diaveg_coeff = np.append(diaveg_coeff, diaveg_calc)
                stemheight_coeff = np.append(stemheight_coeff, stemheight_calc)
                
        
    return drag_coeff, rnveg_coeff, diaveg_coeff, stemheight_coeff

def newCalcDraginLoopCdveg(model_dfm, xzyz_cell_number, index_veg_cel,list_read_subset):
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
    
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen
        
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    
    drag_coeff = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))

    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 0:
            Cd_calc_bare = 0.005
            drag_coeff = np.append(drag_coeff, Cd_calc_bare)

        else:
            # calculate cell_area and water_depth
            cell_area = ba[row] # sudah sesuai karena boundary ada di array paling akhir
            water_depth = hs[row]# sama seperti di atas
            # sometimes there is a condition where the cell is dry
            # therefore, if it is dry I set the drag as Cd_no
            if water_depth == 0:
                Cd_calc_bare = 0.005
                drag_coeff = np.append(drag_coeff, Cd_calc_bare)

            else:
                # print(row)
                 
                read_data_subset = list_read_subset[row]
                
    

                # Volume of the aerial roots
                N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-read_data_subset['dbh_cm']))))
                          
                # if function to adjust the h_pneu for Avicennia
                # this equation is from Du, Qin, et al. 2021
                if water_depth < h_pneu:
                    Vm_pneu = np.pi*d_pneu*water_depth/12
                else:
                    Vm_pneu = np.pi*d_pneu*h_pneu/12
            
                Vm_pneu_total = Vm_pneu*N_pneu #m^3
                Am_pneu_total = d_pneu*h_pneu*N_pneu #m^2
            
                # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
                # height is adjusted from the water depth
                Vm_trunk = np.pi/4*(read_data_subset['dbh_cm']/100)**2*water_depth #m^3 
                Am_trunk = (read_data_subset['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am = Am_pneu_total + Am_trunk#m^2
            
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V-Vm)/Am
                # Cd_no = 0.005 # assume drag coefficient without the presence of vegetation
                e_drag = 5 # set to 5m to obtain realistic value for Cd
            
                Cd_calc = (Cd_no + (e_drag/L)).sum()
                # TODO this is just for reminder of the commented script
                # is it really needed to calculate the weighted drag?
                # Cd_calc_weighted = Cd_calc*(Vm.sum()/V) + Cd_no * ((V-Vm.sum())/V)
                
                # append the calculated cd_veg to the drag_coeff list
                # drag_coeff.append(cd_veg)
                drag_coeff = np.append(drag_coeff, Cd_calc)                
        
    return drag_coeff

def New_clipSHPcreateXLSfromGPD(file_tile, save_tiled_trees, shp_source, species_name):
    """Clip the Master Trees and Create XLS from the GPD read
    a modified function from clipSHPcreateXLSfromGPD
    
    By delete the rbh_m calculation, since the rbh_m has been calculated by
    MesoFON
          
    Parameters
    ---------
    file_tile = path to the tile for looping
    
    save_tiled_trees = path to save the tiled trees and xls
    
    shp_source =  master trees
    
    species_name = species name
    
    a0, b0, a137, b137 = Avicennia Marina Parameters

    Returns
    ---------
    XLS files for MesoFON run
    
    """

    for filepath in glob.iglob(file_tile):
        # calculate the clipping with geopandas
        # gp_point = shp_source
        gp_point= gpd.read_file(shp_source)
        # your_clip = os.path.join(r"D:\Git\d3d_meso\Model-Exchange\Initialization\tile_0_20.shp")
        gp_clip= gpd.read_file(filepath)
    
        tree_point = gp_point.clip(gp_clip)
        # tree_point.to_file(save_loc)
        if tree_point.size > 0:
            tree_point.to_file(os.path.join(save_tiled_trees, Path(filepath).stem + '.shp')) # save tiled shp trees
        # create the XLS file based on the clipped shp
        posX = tree_point['coord_x']
        posY = tree_point['coord_y']
        height_m = tree_point['height_m'] #change for independent files or standardized the shp
        rbh_m = tree_point['rbh_m']
        id_id = np.arange(len(posX))
        speciesName = np.ones(len(posX))*1 #if only one species is recorded, however it is better to place this in shapefile
        types_species = id_id+1 # this variable starts from 1 to N+1
        # shiftedBelowPos = np.ones(len(posX))*1
        age = tree_point['age']
            
        # Create Panda Dataframe with Header
        df = pd.DataFrame({'id': id_id,
                           'speciesName' : speciesName,
                           'type' : types_species,
                           'posX' : posX,
                           'posY' : posY,
                           'shiftedPosX' : posX,
                           'shiftedPosY' : posY,
                           'shiftedBelowPosX' : speciesName,
                           'shiftedBelowPosY' : speciesName,
                           'age' : age,
                           'rbh_m' : rbh_m,
                           'newRbh_m' : rbh_m,
                           'height_m' : height_m,
                           'newHeight_m' : height_m})
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        xls_name = Path(filepath).stem+'_trees_input'+'.xls'
        xls_loc = os.path.join(save_tiled_trees , xls_name)
        
        # Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
        # the default header and skip one row to allow us to insert a user defined
        # header.
        startrowis = 4
        
        with pd.ExcelWriter(xls_loc) as writer:
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=startrowis, header=True)
        
        rb = open_workbook(xls_loc)
        wb = copy(rb)
        # Insert the value in worksheet
        s = wb.get_sheet(0)
        s.write(0,0,'def')
        s.write(1,0, 1)
        s.write(1,1, species_name)
        s.write(2,0, '/def')
        s.write(3,0, 'values')
            # position_last = start the data + 1 (since it contains header) 
        #  + length of data
        position_last = startrowis + 1 + len(posX)
        s.write(position_last,0, '/values')
        wb.save(xls_loc)

def SalNew_func_createRaster4MesoFON(concave_path,save_tiled_env, filename, dir_out, end_name):
    dir_coupling = concave_path
    
    aa = glob.iglob(os.path.join(dir_coupling, filename+'.tif'))
    bb = glob.iglob(os.path.join(dir_out, 'tile_*.tif'))
        
    for (filepatd, filepatc) in itertools.zip_longest(aa,bb):
        # copy and rename the new raster for surv
        raster_surv = filepatd
        raster_surv_name = Path(filepatc).stem+end_name
        target_ras_surv0 = os.path.join(save_tiled_env,raster_surv_name+'0.tif')
        target_ras_surv1 = os.path.join(save_tiled_env,raster_surv_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_surv, os.path.join(save_tiled_env,raster_surv_name+'.tif'))
        shutil.copyfile(raster_surv, target_ras_surv0)
        shutil.copyfile(raster_surv, target_ras_surv1)

def SalNew_Sal_func_createRaster4MesoFON(concave_path,save_tiled_env, filename, dir_out, end_name, val_no_data_sal):
    dir_coupling = concave_path
    
    def changeNoData(datatif,value):
        maskfile = gdal.Open(datatif, gdalconst.GA_Update)
        maskraster = maskfile.ReadAsArray()
        maskraster = np.where((maskraster > -999), maskraster, value ) # yg lebih dari -999 diganti jadi 0
        maskband = maskfile.GetRasterBand(1)
        maskband.WriteArray( maskraster )
        maskband.FlushCache()
        
    aa = glob.iglob(os.path.join(dir_coupling, filename+'.tif'))
    bb = glob.iglob(os.path.join(dir_out, 'tile_*.tif'))
        
    for (filepatd, filepatc) in itertools.zip_longest(aa,bb):
        # copy and rename the new raster for surv
        raster_surv = filepatd
        changeNoData(raster_surv,val_no_data_sal)
        raster_surv_name = Path(filepatc).stem+end_name
        target_ras_surv0 = os.path.join(save_tiled_env,raster_surv_name+'0.tif')
        target_ras_surv1 = os.path.join(save_tiled_env,raster_surv_name+'1.tif')
        # do the copy-paste
        shutil.copyfile(raster_surv, os.path.join(save_tiled_env,raster_surv_name+'.tif'))
        shutil.copyfile(raster_surv, target_ras_surv0)
        shutil.copyfile(raster_surv, target_ras_surv1)

def calcLevelCell (xyzw_cell_number, model_dfm, index_veg_cel, xyzw_nodes, xk, yk, read_data):
    bulk_density = 1.1 ##gr/cm^3
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    
    bl_val = np.empty((model_dfm.get_var('bl').shape[0],0))
    for row in range(len(xyzw_cell_number)):
        if index_veg_cel[row] == 1:
            cell_area = ba[row]
            
            position = xyzw_cell_number[row,2].astype(int)
                            
            nodes_data = ma.compressed(xyzw_nodes[position][xyzw_nodes[position].mask == False]).astype(int)# select only the valid data (unmasked / false)
            nodes_pos = np.block([[xk[nodes_data-1]],[yk[nodes_data-1]]]) # substracted to 1 in order to adjust the 0-based position in python
            # Find the min max of each x,y coordinate
            # create the list of x_min-x_max and y_min-y_max
            x_range = [np.min(nodes_pos[0]), np.max(nodes_pos[0])]
            y_range = [np.min(nodes_pos[1]), np.max(nodes_pos[1])]
                            
            # subsetting pandas 
            read_data_subset = read_data[(read_data['GeoRefPosX'] >= x_range[0]) & 
                                         (read_data['GeoRefPosX'] <= x_range[1])]
            read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= y_range[0]) & 
                                                (read_data_subset['GeoRefPosY'] <= y_range[1])]
            
            bgb = 1.28 *  read_data_subset['dbh_cm']**1.17
            bgb_kg_total = bgb.sum()
            bgb_vol = bgb_kg_total * 1000 / bulk_density * 1E-6# m^3
            bgb_level_veg = bgb_vol / cell_area # normalized bed level increase 
            
            bl_val = np.append(bl_val, bgb_level_veg)
        else:
            bgb_vol_bare = 0
            bl_val = np.append(bl_val, bgb_vol_bare)
    
    return bl_val

def calcLevelCellCdveg (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data):
    bulk_density = 1.1 ##gr/cm^3
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    data_verts = data_verts = ugrid_all.verts
    
    bl_val = np.empty((model_dfm.get_var('bl').shape[0],0))
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 1:
            cell_area = ba[row]
            
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                         (read_data['GeoRefPosX'] <= max(aa[:,0]))]
            read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                                (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
            
            bgb = 1.28 *  read_data_subset['dbh_cm']**1.17
            bgb_kg_total = bgb.sum()
            bgb_vol = bgb_kg_total * 1000 / bulk_density * 1E-6# m^3
            bgb_level_veg = bgb_vol / cell_area # normalized bed level increase 
            
            bl_val = np.append(bl_val, bgb_level_veg)
        else:
            bgb_vol_bare = 0
            bl_val = np.append(bl_val, bgb_vol_bare)
    
    return bl_val

def definePropVeg(xzyz_cell_number, model_dfm, ugrid_all, index_veg_cel, read_data, addition_veg):
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    ## initialisation of vegetation variables

    data_verts = ugrid_all.verts
    
    rnveg_coeff = np.zeros((model_dfm.get_var('Cdvegsp').shape[0],0))
    diaveg_coeff = np.zeros((model_dfm.get_var('Cdvegsp').shape[0],0))
    stemheight_coeff = np.zeros((model_dfm.get_var('Cdvegsp').shape[0],0))
    
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 0:
            rnveg_coeff = np.append(rnveg_coeff, 0)
            diaveg_coeff = np.append(diaveg_coeff, 0)
            stemheight_coeff = np.append(stemheight_coeff, 0)
        else:
            # calculate cell_area and water_depth
            cell_area = ba[row] # sudah sesuai karena boundary ada di array paling akhir

            # print(row)
            # find the position based on the cell number
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                         (read_data['GeoRefPosX'] <= max(aa[:,0]))]
            read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                                (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
           
            rnveg_calc = read_data_subset.Height_cm.notnull().sum()/cell_area
            diaveg_calc = read_data_subset.dbh_cm.quantile(0.8)/100
            stemheight_calc = read_data_subset.Height_cm.quantile(0.8)/100
            
            rnveg_coeff = np.append(rnveg_coeff, rnveg_calc)
            diaveg_coeff = np.append(diaveg_coeff, diaveg_calc)
            stemheight_coeff = np.append(stemheight_coeff, stemheight_calc)
     
    rnveg_coeff = np.append(rnveg_coeff, addition_veg)
    diaveg_coeff = np.append(diaveg_coeff, addition_veg)
    stemheight_coeff = np.append(stemheight_coeff, addition_veg)
    
    return rnveg_coeff, diaveg_coeff, stemheight_coeff