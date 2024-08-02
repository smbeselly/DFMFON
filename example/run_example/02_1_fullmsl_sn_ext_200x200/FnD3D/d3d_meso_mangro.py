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
from scipy.integrate import quad
from scipy.interpolate import interp2d

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
    
    for row in np.arange(len(data_verts[:,0,0])):
        aa = data_verts[row,:,:]
                
        ab_subset = df_xzyz.loc[(((df_xzyz[0] >= min(aa[:,0])) 
                            & (df_xzyz[0] <= max(aa[:,0]))) 
                                          & ((df_xzyz[1] >= min(aa[:,1])) 
                                              & (df_xzyz[1] <= max(aa[:,1]))))]
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
    read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= x_range[0]) 
                            & (read_data['GeoRefPosX'] <= x_range[1])) 
                                          & ((read_data['GeoRefPosY'] >= y_range[0]) 
                                             & (read_data['GeoRefPosY'] <= y_range[1])))]

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
    
    if WoOfin.size == 0:
        fromheight = [0,0]
        Pvalue = [0,0]
        print('WoOfin size == 0: WoO Zero')
    elif max(step) == 0:
        fromheight = [0,0]
        Pvalue = [0,0]
        print('Step == 0: WoO Zero')
    else:
        ### probability of survival for different elevations, P success ####
        d1 = WoOfin[(WoOfin[:,0] == 1)]
        step = np.arange(np.ceil(np.sqrt(len(d1)))+1)
        ## sqrt of N,round up, define number of steps for probability analysis
        ## along height gradient
        
        Pvalue = np.zeros(int(max(step)))
        fromheight = np.zeros(int(max(step)))
        toheight = np.zeros(int(max(step)))
        PvalueNo = np.zeros(int(max(step)))
        startel = min(WoOfin[:,1])
        stepsize = (max(WoOfin[:,1]) - min(WoOfin[:,1]))/max(step)

        ## calculate size of steps probability calculation
        for i in range(int(max(step))):
            ii = i+1
            Pselect2 = 0
            Pselect1 = 0
            Pselect2 = WoOfin[(WoOfin[:,1] >= (startel + (stepsize*ii)))]
            Pselect1 = Pselect2[(Pselect2[:,1] <= (startel + (stepsize*(1+ii))))]  
            if Pselect1.size > 0:
                Pvalue[i] = sum(Pselect1[:,0])/len(Pselect1[:,0])
                PvalueNo[i] = (sum(Pselect1[:,0])/len(Pselect1[:,0])) * sum(Pselect1[:,0])
            else:
                Pvalue[i] = 0
                PvalueNo[i] = 0
            # Pvalue[i] = sum(Pselect1[:,0])/len(Pselect1[:,0])
            # PvalueNo[i] = (sum(Pselect1[:,0])/len(Pselect1[:,0])) * sum(Pselect1[:,0])
            fromheight[i] = startel + (stepsize * (ii-1))
            toheight[i] = startel + (stepsize * (ii))
    
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

def n_calcAgeCoupling0(read_data, mature_sapl):
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
    master_trees_np = mature_sapl[['GeoRefPosX','GeoRefPosY']].to_numpy()
    # assign the correct cell number to the matrix (cell number is python based)
    # if the cell number is read as 0, the real cell number in netCDF is 1
    #source: https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
    age_read_data = np.empty((len(read_data),0))
    for row in range(len(read_data)): #source: https://www.pluralsight.com/guides/numpy-arrays-iterating
        pt = [read_data_np[row,:]]
        index = spatial.KDTree(master_trees_np).query(pt)
        age_read_data = np.append(age_read_data,mature_sapl['Age'][int(index[1])])
        # age_read_data.append(master_trees['age'][int(index[1])])

    age_coupling0 = read_data['tick']+age_read_data
    age_coupling0 = age_coupling0.reset_index(drop=True)
    # age_coupling0 = age_coupling.fillna(0.25)
    
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
    
def createPointSHPmulti(read_data, age_coupling, path_shp, EPSG_Project):
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
                  read_data['Height_cm']/100,age_coupling,read_data['dbh_cm']/200,
                  read_data['Species']]
    data_point_header = ['coord_x','coord_y','height_m','age','rbh_m','Species']
    data_point_df = pd.concat(data_point, axis=1, keys=data_point_header)
    # define schema
    schema = {
        'geometry':'Point',
        'properties':{'coord_x':'float',
                      'coord_y':'float',
                      'height_m':'float',
                      'age':'float',
                      'rbh_m':'float',
                      'Species':'str'}
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
                           'rbh_m' : row.rbh_m,
                           'Species': row.Species},
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
        read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= x_range[0]) 
                            & (read_data['GeoRefPosX'] <= x_range[1])) 
                                          & ((read_data['GeoRefPosY'] >= y_range[0]) 
                                             & (read_data['GeoRefPosY'] <= y_range[1])))]
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
                try:
                    prep_del = prep_out
                    del(prep_out)
                    os.remove(str(prep_del))
                except:
                    pass
                
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
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= x_range[0]) 
                            & (read_data['GeoRefPosX'] <= x_range[1])) 
                                          & ((read_data['GeoRefPosY'] >= y_range[0]) 
                                             & (read_data['GeoRefPosY'] <= y_range[1])))]
        
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
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
        
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
                read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= x_range[0]) 
                            & (read_data['GeoRefPosX'] <= x_range[1])) 
                                          & ((read_data['GeoRefPosY'] >= y_range[0]) 
                                             & (read_data['GeoRefPosY'] <= y_range[1])))]
               
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
    
    def calc_ab(HR, x2):
        if HR > 2.24:
            b = 0
            a = -(HR/x2**2)
        else:
            b = np.tan((20.88*HR - 46.87)*np.pi/180)
            a = -((b*x2+HR)/x2**2)
        
        return a, b

    def calc_F_vol(a,b,x):
        F = (((2*a*x + b) * np.sqrt(1 + (2*a*x + b)**2)) / 4*a) + \
            ((np.log((2*a*x + b) + np.sqrt(1 + (2*a*x + b)**2)))/4*a)
        
        return F
    
    def calc_L_vol(x2, HR, a, b, D):
        if D < HR:
            x1 = (-b-np.sqrt(b**2-4*a*(HR-D)))/2*a
            F1 = calc_F_vol(a,b,x1)
            F2 = calc_F_vol(a,b,x2)
            L = F2-F1
        else:
            F2 = calc_F_vol(a,b,x2)
            F0 = calc_F_vol(a,b,0)
            L = F2-F0
        
        return L
    
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
                read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
                
                avicennia = read_data_subset.loc[read_data_subset['Species']=='Avicennia_marina'].copy() # for avicennia
                RhDFrm = read_data_subset.loc[read_data_subset['Species']=='Rhizopora_apiculata'].copy() # for Rhizopora
                
                ## Calculate Avicennia Contribution
                # Volume of the aerial roots
                N_pneu = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-avicennia['dbh_cm']))))
                          
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
                Vm_trunk = np.pi/4*(avicennia['dbh_cm']/100)**2*water_depth #m^3 
                Am_trunk = (avicennia['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm_avc = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am_avc = Am_pneu_total + Am_trunk#m^2
                
                ## Contribution of Rhizopora; Ohira, et.al (2013)
                
                RhDFrm['HRmax'] = (7.56*RhDFrm['dbh_cm']/100) + 0.5
                RhDFrm['N'] = np.ceil(3.15*RhDFrm['HRmax']**2 + 5.3*RhDFrm['HRmax'] + 0.114)
                RhDFrm['intervalRoot']= (RhDFrm['HRmax']-0.25) / (RhDFrm['N'] -1) # this is the interval of each root shoots + 0.25
                RhDFrm['HR'] = 0.25+RhDFrm['intervalRoot']
                RhDFrm['x2'] = 1.88*RhDFrm['HR'] - 0.432
                RhDFrm['rootDia'] = 0.04*RhDFrm['dbh_cm']/100 + 0.005*RhDFrm['HR'] + 0.024
                RhDFrm['thetaRoot'] = 20.88*RhDFrm['HR'] - 46.87 # in degrees 
                
                # HRmax = RhDFrm['HRmax'].to_numpy()
                HR = RhDFrm['HR'].to_numpy()
                x2 = RhDFrm['x2'].to_numpy()
                rootDia = RhDFrm['rootDia'].to_numpy()
                RhDBH = RhDFrm['dbh_cm'].to_numpy()/100
                
                a = []
                b = []
                for i in np.arange(len(HR)):
                    acalc, bcalc = calc_ab(HR[i], x2[i])
                    
                    a = np.append(a, acalc)
                    b = np.append(b, bcalc)
                
                L = []  ## Ini yang butuh water level
                for i in np.arange(len(x2)):
                    Lcalc = calc_L_vol(x2[i], HR[i], a[i], b[i], water_depth)
                    
                    L = np.append(L, Lcalc)
                
                Vol_Root = L*(rootDia/2)**2*np.pi
                Area_Root = L*rootDia
                Vm_stem = []
                Am_stem = []
                for i in np.arange(len(RhDBH)):
                    if RhDBH[i] > water_depth:
                        Vm_stemp = 0
                        Am_stemp = 0
                    else:
                        Vm_stemp = np.pi/4*(RhDBH[i])**2*water_depth #m^3 
                        Am_stemp = RhDBH[i]*water_depth #m^3 
                    
                    Vm_stem = np.append(Vm_stem, Vm_stemp)
                    Am_stem = np.append(Am_stem, Am_stemp)
                
                # sum of the obstacle volume
                Vm_rhiz = pd.Series(Vol_Root + Vm_stem) #m^3
                # sum of the obstacle area
                Am_rhiz = pd.Series(Area_Root + Am_stem)#m^2
                
                Vm_tot = Vm_avc.append(Vm_rhiz, ignore_index=True)
                Am_tot = Am_avc.append(Am_rhiz, ignore_index=True)
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V - Vm_tot) / Am_tot
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
    
    def calc_ab(HR, x2):
        if HR > 2.24:
            b = 0
            a = -(HR/x2**2)
        else:
            b = np.tan((20.88*HR - 46.87)*np.pi/180)
            a = -((b*x2+HR)/x2**2)
        
        return a, b

    def calc_F_vol(a,b,x):
        F = (((2*a*x + b) * np.sqrt(1 + (2*a*x + b)**2)) / 4*a) + \
            ((np.log((2*a*x + b) + np.sqrt(1 + (2*a*x + b)**2)))/4*a)
        
        return F
    
    def calc_L_vol(x2, HR, a, b, D):
        if D < HR:
            x1 = (-b-np.sqrt(b**2-4*a*(HR-D)))/2*a
            F1 = calc_F_vol(a,b,x1)
            F2 = calc_F_vol(a,b,x2)
            L = F2-F1
        else:
            F2 = calc_F_vol(a,b,x2)
            F0 = calc_F_vol(a,b,0)
            L = F2-F0
        
        return L
    
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
                
                avicennia = read_data_subset.loc[read_data_subset['Species']=='Avicennia_marina'].copy() # for avicennia
                RhDFrm = read_data_subset.loc[read_data_subset['Species']=='Rhizopora_apiculata'].copy() # for Rhizopora

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
                Am_trunk = (avicennia['dbh_cm']/100)*water_depth #m^3 
            
                # sum of the obstacle volume
                Vm_avc = Vm_pneu_total + Vm_trunk #m^3
                # sum of the obstacle area
                Am_avc = Am_pneu_total + Am_trunk#m^2
                
                ## Contribution of Rhizopora
                
                RhDFrm['HRmax'] = (7.56*RhDFrm['dbh_cm']/100) + 0.5
                RhDFrm['N'] = np.ceil(3.15*RhDFrm['HRmax']**2 + 5.3*RhDFrm['HRmax'] + 0.114)
                RhDFrm['intervalRoot']= (RhDFrm['HRmax']-0.25) / (RhDFrm['N'] -1) # this is the interval of each root shoots + 0.25
                RhDFrm['HR'] = 0.25+RhDFrm['intervalRoot']
                RhDFrm['x2'] = 1.88*RhDFrm['HR'] - 0.432
                RhDFrm['rootDia'] = 0.04*RhDFrm['dbh_cm']/100 + 0.005*RhDFrm['HR'] + 0.024
                RhDFrm['thetaRoot'] = 20.88*RhDFrm['HR'] - 46.87 # in degrees 
                
                # HRmax = RhDFrm['HRmax'].to_numpy()
                HR = RhDFrm['HR'].to_numpy()
                x2 = RhDFrm['x2'].to_numpy()
                rootDia = RhDFrm['rootDia'].to_numpy()
                RhDBH = RhDFrm['dbh_cm'].to_numpy()/100
                
                a = []
                b = []
                for i in np.arange(len(HR)):
                    acalc, bcalc = calc_ab(HR[i], x2[i])
                    
                    a = np.append(a, acalc)
                    b = np.append(b, bcalc)
                
                L = []  ## Ini yang butuh water level
                for i in np.arange(len(x2)):
                    Lcalc = calc_L_vol(x2[i], HR[i], a[i], b[i], water_depth)
                    
                    L = np.append(L, Lcalc)
                
                Vol_Root = L*(rootDia/2)**2*np.pi
                Area_Root = L*rootDia
                Vm_stem = []
                Am_stem = []
                for i in np.arange(len(RhDBH)):
                    if RhDBH[i] > water_depth:
                        Vm_stemp = 0
                        Am_stemp = 0
                    else:
                        Vm_stemp = np.pi/4*(RhDBH[i])**2*water_depth #m^3 
                        Am_stemp = RhDBH[i]*water_depth #m^3 
                    
                    Vm_stem = np.append(Vm_stem, Vm_stemp)
                    Am_stem = np.append(Am_stem, Am_stemp)
                
                # sum of the obstacle volume
                Vm_rhiz = pd.Series(Vol_Root + Vm_stem) #m^3
                # sum of the obstacle area
                Am_rhiz = pd.Series(Area_Root + Am_stem)#m^2
                
                Vm_tot = Vm_avc.append(Vm_rhiz, ignore_index=True)
                Am_tot = Am_avc.append(Am_rhiz, ignore_index=True)
                # Volume of the control water
            
                V = cell_area*water_depth
            
                # Characteristic Length (m)
                L = (V - Vm_tot) / Am_tot
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
        speciesName = tree_point['Species']
        types_species = id_id+1 # this variable starts from 1 to N+1
        shiftedBelowPos = np.ones(len(posX))*1
        age = tree_point['age']
            
        # Create Panda Dataframe with Header
        df = pd.DataFrame({'id': id_id,
                           'speciesName' : speciesName,
                           'type' : types_species,
                           'posX' : posX,
                           'posY' : posY,
                           'shiftedPosX' : posX,
                           'shiftedPosY' : posY,
                           'shiftedBelowPosX' : shiftedBelowPos,
                           'shiftedBelowPosY' : shiftedBelowPos,
                           'age' : age,
                           'rbh_m' : rbh_m,
                           'newRbh_m' : rbh_m,
                           'height_m' : height_m,
                           'newHeight_m' : height_m})
        for nn in range(len(species_name)):
            df.loc[df["speciesName"] == species_name[nn], "speciesName"] = nn+1
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        xls_name = Path(filepath).stem+'_trees_input'+'.xls'
        xls_loc = os.path.join(save_tiled_trees , xls_name)
        
        # Convert the dataframe to an XlsxWriter Excel object. Note that we turn off
        # the default header and skip one row to allow us to insert a user defined
        # header.
        startrowis = 3+len(species_name)
    
        with pd.ExcelWriter(xls_loc) as writer:
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=startrowis, header=True)
        
        rb = open_workbook(xls_loc)
        wb = copy(rb)
        # Insert the value in worksheet
        s = wb.get_sheet(0)
        s.write(0,0,'def')
            
        for nn in range(len(species_name)):
            s.write(nn+1,0, nn+1)
            s.write(nn+1,1, species_name[nn])
        s.write(len(species_name)+1,0, '/def')
        s.write(len(species_name)+2,0, 'values')
            
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
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= x_range[0]) 
                            & (read_data['GeoRefPosX'] <= x_range[1])) 
                                          & ((read_data['GeoRefPosY'] >= y_range[0]) 
                                             & (read_data['GeoRefPosY'] <= y_range[1])))]
            
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
    data_verts = ugrid_all.verts
    
    bl_val = np.empty((model_dfm.get_var('bl').shape[0],0))
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 1:
            cell_area = ba[row]
            
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
            
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
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
           
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

def n_definePropVeg(avicennia_pr, RhDFrm_pr, vol_func, model_dfm):
    avicennia = avicennia_pr.copy()
    RhDFrm = RhDFrm_pr.copy()
    
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    
    lookup_area = pd.DataFrame({'cell_number':list(range(len(ba))),
                               'cell_area':ba,
                               })
    try:
        # avicennia
        val_ba_avc = np.vstack(avicennia['cell_number'].values)
        bound_ba_avc = lookup_area['cell_number'].values == val_ba_avc
        avicennia['cell_area'] = np.dot(bound_ba_avc, lookup_area['cell_area'])
    except:
        print('no avicennia')
    try:
        # rhizopora
        val_ba_rzp = np.vstack(RhDFrm['cell_number'].values)
        bound_ba_rzp = lookup_area['cell_number'].values == val_ba_rzp
        RhDFrm['cell_area'] = np.dot(bound_ba_rzp, lookup_area['cell_area'])
    except:
        print('no avicennia')
        
    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    
    # define HRmax
    RhDFrm['HRmax'] = ((7.56*RhDFrm['dbh_cm']/100) + 0.5).astype(float)
    
    ## define water depth as the limit
    avicennia['water_depth'] = h_pneu
    RhDFrm['water_depth'] = RhDFrm['HRmax']
    
    ## calculate contribution from avicennia
    try:
        cond_avc = [
            (avicennia['water_depth'].values < h_pneu),
            (avicennia['water_depth'].values >= h_pneu)
            ]
        val_avc = [
            np.pi*d_pneu*avicennia['water_depth'].values/12,
            np.pi*d_pneu*h_pneu/12
            ]
    
        avicennia['Vm_pneu'] = np.select(cond_avc, val_avc, default='NA').astype(float)
    
        avicennia['Vm_pneu_total'] = avicennia['Vm_pneu']*avicennia['N_pneu']
        # avicennia['Am_pneu_total'] = d_pneu*h_pneu*avicennia['N_pneu']
    
        # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
        # height is adjusted from the water depth
        avicennia['Vm_trunk'] = np.pi/4*(avicennia['dbh_cm']/100)**2*avicennia['water_depth']#m^3 
        # avicennia['Am_trunk'] = (avicennia['dbh_cm']/100)*avicennia['water_depth'] #m^3 
    
        # sum of the obstacle volume
        avicennia['Vm_avc'] = avicennia['Vm_pneu_total'] + avicennia['Vm_trunk'] #m^3
        # sum of the obstacle area
        # avicennia['Am_avc'] = avicennia['Am_pneu_total'] + avicennia['Am_trunk']#m^2  
        avicennia['nmang_eq_avc'] = np.round(avicennia['Vm_avc']  / avicennia['Vm_trunk']) + 1
    except:
        print('no avicennia')
        
    ## calculate contribution from rhizopora
    try:           
        # RhDFrm['Vm_rhiz'] = g_func(RhDFrm['water_depth'], RhDFrm['dbh_cm']) #water depth(m), dbh(cm)
        
        data_array = RhDFrm[["water_depth", "dbh_cm"]].to_numpy()
        def calc_vol(val):
            return  vol_func(val[0], val[1])
        # def calc_are(val):
        #     return  area_func(val[0], val[1])
        
        RhDFrm['Vm_rhiz'] = np.apply_along_axis(calc_vol, axis=1, arr=data_array)
        RhDFrm['Vm_trunk'] = np.pi/4*(RhDFrm['dbh_cm']/100)**2*RhDFrm['water_depth']#m^3 
        RhDFrm['stilt_roots'] = RhDFrm['Vm_rhiz'] - RhDFrm['Vm_trunk']
        RhDFrm['nmang_eq_rzp'] = np.round(RhDFrm['stilt_roots'] / RhDFrm['Vm_trunk']) + 1
        # RhDFrm['Am_rhiz'] = np.apply_along_axis(calc_are, axis=1, arr=data_array)
    except:
        print('no rhizopora')
        
    try:
        avicennia_sum = avicennia.groupby(by='cell_number')['nmang_eq_avc'].sum()
        avicennia_dia = avicennia.groupby(by='cell_number')['dbh_cm'].quantile(q=0.8)
        avicennia_height = avicennia.groupby(by='cell_number')['Height_cm'].quantile(q=0.8)
    except:
        print('no avicennia')
    try:
        rhizopora_sum = RhDFrm.groupby(by='cell_number')['nmang_eq_rzp'].sum()
        rhizopora_dia = RhDFrm.groupby(by='cell_number')['dbh_cm'].quantile(q=0.8)
        rhizopora_height = RhDFrm.groupby(by='cell_number')['Height_cm'].quantile(q=0.8)
        
        veg_contribute = pd.concat([avicennia_sum, rhizopora_sum], axis=1)
        
    except:
        print('no rhizopora')
    
    veg_prop = pd.concat([lookup_area, veg_contribute], axis=1)
    veg_prop['rnveg_calc'] = (veg_prop['nmang_eq_avc'] + veg_prop['nmang_eq_rzp']) / veg_prop['cell_area']
    rnveg_calc = veg_prop['rnveg_calc'].fillna(0).to_numpy()
    
    diaveg_is = pd.concat([avicennia_dia, rhizopora_dia], axis=1)
    diaveg_is.columns = ['avicennia_dbh_cm', 'rhizopora_dbh_cm']
    diaveg_precalc = diaveg_is.quantile(q=0.8, axis=1)
    
    stemheight_is = pd.concat([avicennia_height, rhizopora_height], axis=1)
    stemheight_is.columns = ['avicennia_height_cm', 'rhizopora_height_cm']
    stemheight_precalc = stemheight_is.quantile(q=0.8, axis=1)
    
    diaveg_calc = pd.concat([lookup_area, diaveg_precalc], axis=1)
    diaveg_calc = diaveg_calc.fillna(0)
    diaveg_calc.columns = ['cell_number', 'cell_area', 'diaveg']
    diaveg_calc_cell = (diaveg_calc['diaveg']/100).to_numpy()
    
    stemheight_calc = pd.concat([lookup_area, stemheight_precalc], axis=1)
    stemheight_calc = stemheight_calc.fillna(0)
    stemheight_calc.columns = ['cell_number', 'cell_area', 'stemheight']
    stemheight_calc_cell = (stemheight_calc['stemheight']/100).to_numpy()
    
    return rnveg_calc, diaveg_calc_cell, stemheight_calc_cell

def Calc_WoO(water_level, model_dfm, MorFac, coupling_period, woo_inun, LLWL):
    import copy           
    ### Calculate the WoO value
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

    ### Calculate the probability based on the WoO and store the value in each cell
    # find median value of the h_wl for each cell number
    med_h_wl = np.median(h_wl, axis=1) 

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
            
    return med_h_wl, surv_val

def fill_elev_base (model_dfm, ugrid_all, xzyz_cell_number, index_veg_cel, read_data):

    data_verts = ugrid_all.verts
    
    # bl_val = np.empty((model_dfm.get_var('bl').shape[0],0))
    bl_val = model_dfm.get_var('bl')
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 1: 
            # print(row)
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
            
            # read_data_subset['elev_base'] = bl_val[row]
            pos = read_data_subset.index.values.tolist()
            read_data.loc[pos,'elev_base'] = bl_val[row]
        else:
            pass

def fill_elev_sdlg (bed_lvl_init, ugrid_all, xzyz_cell_number, index_veg_cel, read_data):

    data_verts = ugrid_all.verts
    
    # bl_val = np.empty((model_dfm.get_var('bl').shape[0],0))
    for row in range(len(xzyz_cell_number)):
        if index_veg_cel[row] == 1: 
            # print(row)
            position = xzyz_cell_number[row,2].astype(int)
            aa = data_verts[position,:,:]
            # subsetting pandas 
            read_data_subset = read_data.loc[(((read_data['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data['GeoRefPosX'] <= max(aa[:,0]))) 
                                          & ((read_data['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data['GeoRefPosY'] <= max(aa[:,1]))))]
            
            # read_data_subset['elev_base'] = bl_val[row]
            pos = read_data_subset.index.values.tolist()
            read_data.loc[pos,'elev_base'] = bed_lvl_init[row]
        else:
            pass

## pakai function ini di atas sendiri, muncul di setiap coupling saja
def add_cell_number(ugrid_all, xzyz_cell_number, read_data):
    # this function is to predefine the subset before calculate drag
    read_data_copy = read_data.copy()
    read_data_copy['cell_number'] = np.nan
    data_verts = ugrid_all.verts    
    for row in range(len(xzyz_cell_number)):
        position = xzyz_cell_number[row,2].astype(int)
        aa = data_verts[position,:,:]
        
        read_data_copy['cell_id'] = ((read_data_copy['GeoRefPosX'] >= min(aa[:,0])) 
                            & (read_data_copy['GeoRefPosX'] <= max(aa[:,0]))) & ((read_data_copy['GeoRefPosY'] >= min(aa[:,1])) 
                                             & (read_data_copy['GeoRefPosY'] <= max(aa[:,1])))
        read_data_copy['cell_number'] = np.where(read_data_copy['cell_id'] == True, row, read_data_copy['cell_number'])
    
    return read_data_copy

def newCalcDraginVectCdveg(model_dfm, ugrid_all, xzyz_cell_number, read_data, index_veg_cel):
       
    # def calc_F_vol(a,b,x):
    #         F = (((2*a*x + b) * np.sqrt(1 + (2*a*x + b)**2)) / 4*a) + \
    #             ((np.log((2*a*x + b) + np.sqrt(1 + (2*a*x + b)**2)))/4*a)
            
    #         return F
        
    # def calc_L_vol(x2, HR, a, b, D):
    #     if D < HR:
    #         x1 = (-b-np.sqrt(b**2-4*a*(HR-D)))/2*a
    #         F1 = calc_F_vol(a,b,x1)
    #         F2 = calc_F_vol(a,b,x2)
    #         L = F2-F1
    #     else:
    #         F2 = calc_F_vol(a,b,x2)
    #         F0 = calc_F_vol(a,b,0)
    #         L = F2-F0
        
    #     return L
    # read_data_copy = read_data.copy()
                  
    read_data_copy = read_data.copy()
    ### define cell area and water depth
    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}
    
    lookup_area = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'cell_area':ba,
                           })
    lookup_hs = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'water_depth':hs,
                           })
    
    val_ba = np.vstack(read_data_copy['cell_number'].values)
    bound_ba = lookup_area['cell_number'].values == val_ba
    read_data_copy['cell_area'] = np.dot(bound_ba, lookup_area['cell_area'])
    
    bound_hs = lookup_hs['cell_number'].values == val_ba
    read_data_copy['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
    
    ### calculate for each species
    ## Avicennia
    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation

    d_pneu = 0.01 # m
    h_pneu = 0.15 # m
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen
    
    avicennia = (read_data_copy.loc[read_data_copy['Species']=='Avicennia_marina']).copy() # for avicennia
    avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    # Volume of the aerial roots
    avicennia['N_pneu'] = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-avicennia['dbh_cm']))))
              
    # if function to adjust the h_pneu for Avicennia
    # this equation is from Du, Qin, et al. 2021
    cond_avc = [
        (avicennia['water_depth'].values < h_pneu),
        (avicennia['water_depth'].values >= h_pneu)
        ]
    val_avc = [
        np.pi*d_pneu*avicennia['water_depth'].values/12,
        np.pi*d_pneu*h_pneu/12
        ]
    
    avicennia['Vm_pneu'] = np.select(cond_avc, val_avc, default='NA').astype(float)
    
    avicennia['Vm_pneu_total'] = avicennia['Vm_pneu']*avicennia['N_pneu']
    avicennia['Am_pneu_total'] = d_pneu*h_pneu*avicennia['N_pneu']
    
    # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
    # height is adjusted from the water depth
    avicennia['Vm_trunk'] = np.pi/4*(avicennia['dbh_cm']/100)**2*avicennia['water_depth']#m^3 
    avicennia['Am_trunk'] = (avicennia['dbh_cm']/100)*avicennia['water_depth'] #m^3 
    
    # sum of the obstacle volume
    avicennia['Vm_avc'] = avicennia['Vm_pneu_total'] + avicennia['Vm_trunk'] #m^3
    # sum of the obstacle area
    avicennia['Am_avc'] = avicennia['Am_pneu_total'] + avicennia['Am_trunk']#m^2  
    
    ### Contribution Rhizopora    
    RhDFrm = (read_data_copy.loc[read_data_copy['Species']=='Rhizopora_apiculata']).copy() # for Rhizopora
    RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)
    
    RhDFrm['HRmax'] = ((7.56*RhDFrm['dbh_cm']/100) + 0.5).astype(float)
    RhDFrm['N'] = np.ceil(3.15*RhDFrm['HRmax']**2 + 5.3*RhDFrm['HRmax'] + 0.114).astype(int)
    RhDFrm['HR_interval'] = (RhDFrm['HRmax'] -0.25)/(RhDFrm['N']-1)
    # create linspace of root height as list
    RhDFrm['HR_i'] = RhDFrm.apply(lambda x: np.linspace(0.25,  x['HRmax'], x['N']), axis=1)
    
    # explode the DataFrame
    RhDFrm_ex = RhDFrm.explode('HR_i')
    RhDFrm_ex['HR_i'] = RhDFrm_ex['HR_i'].astype(float)
    # calculate as normally
    RhDFrm_ex['x2'] = (1.88*RhDFrm_ex['HR_i'])-0.432
    RhDFrm_ex['rootDia'] = (0.04*RhDFrm_ex['dbh_cm']/100) + (0.005*RhDFrm_ex['HR_i']) + 0.024
    RhDFrm_ex['thetaRoot'] = 20.88*RhDFrm_ex['HR_i'] - 46.87 # in degrees
    RhDFrm_ex['b'] = np.where(RhDFrm_ex['HR_i'] > 2.24,
        0,
        np.tan((20.88*RhDFrm_ex['HR_i'] - 46.87)*np.pi/180)
        )
    RhDFrm_ex['a'] = np.where(RhDFrm_ex['HR_i'] > 2.24, 
             -(RhDFrm_ex['HR_i']/RhDFrm_ex['x2']**2), 
             -((RhDFrm_ex['b']*RhDFrm_ex['x2']+RhDFrm_ex['HR_i'])/RhDFrm_ex['x2']**2)
             )
    
    RhDFrm_ex['x1'] = np.where(RhDFrm_ex['water_depth'] < RhDFrm_ex['HR_i'],
        (-RhDFrm_ex['b']-np.sqrt(RhDFrm_ex['b']**2-(4*RhDFrm_ex['a']*(RhDFrm_ex['HR_i']-RhDFrm_ex['water_depth']))))/(2*RhDFrm_ex['a']),
        0
        )
    
    calc_L_array = RhDFrm_ex[["a", "b", "x1", "x2"]].to_numpy()
    def calc_L(x):
        g = lambda f, a,b: np.sqrt(1+((2*a*f)+b)**2)
        return quad(g, x[2],x[3], args=(x[0],x[1]))[0]
    RhDFrm_ex['L'] = np.apply_along_axis(calc_L, axis=1, arr=calc_L_array)
    
    RhDFrm_ex['VR_i'] = RhDFrm_ex['L']*(RhDFrm_ex['rootDia']/2)**2
    RhDFrm_ex['AR_i'] = RhDFrm_ex['L']*RhDFrm_ex['rootDia']
    RhDFrm_grouped = RhDFrm_ex.groupby(by='Unnamed: 0')[['VR_i','AR_i']].sum().to_numpy()
    
    #add the calculated V and A to main DataFrame
    RhDFrm['Vol_Root'] = RhDFrm_grouped[:,0].tolist()
    RhDFrm['Area_Root'] = RhDFrm_grouped[:,1].tolist()
    # cond_rzp = [
    #     (RhDFrm['water_depth'].values < RhDFrm['HR']),
    #     (RhDFrm['water_depth'].values >= RhDFrm['HR'])
    #     ]
    # # check x1
    # try:
    #     answer_x1 = [
    #         (-RhDFrm['b']-np.sqrt(RhDFrm['b']**2-4*RhDFrm['a']*(RhDFrm['HR']-RhDFrm['water_depth'])))/2*RhDFrm['a'],
    #         np.nan
    #         ]
    #     RhDFrm['x1'] = np.select(cond_rzp, answer_x1, default='NA').astype(float)
    # except:
    #     RhDFrm['x1'] = np.nan
    # #check F1
    # try:
    #     answer_F1 = [
    #         calc_F_vol(RhDFrm['a'],RhDFrm['b'],RhDFrm['x1']),
    #         np.nan
    #         ]
        
    #     RhDFrm['F1'] = np.select(cond_rzp, answer_F1, default='NA').astype(float)
    # except:
    #     RhDFrm['F1'] = np.nan
    
    # #check F2
    # RhDFrm['F2'] = calc_F_vol(RhDFrm['a'],RhDFrm['b'],RhDFrm['x2'])
    
    # #check F0
    # try:
    #     answer_F0 = [
    #         np.nan,
    #         calc_F_vol(RhDFrm['a'],RhDFrm['b'],0)
    #         ]
    #     RhDFrm['F0'] = np.select(cond_rzp, answer_F0, default='NA').astype(float)
    # except:
    #     RhDFrm['F0'] = np.nan
    # #Check L
    # answer_L = [
    #     RhDFrm['F2']-RhDFrm['F1'],
    #     RhDFrm['F2']-RhDFrm['F0']
    #     ]
    # RhDFrm['L'] = np.select(cond_rzp, answer_L, default='NA').astype(float)
    
    # RhDFrm['Vol_Root'] = RhDFrm['L']*(RhDFrm['rootDia']/2)**2*np.pi
    # RhDFrm['Area_Root'] = RhDFrm['L']*RhDFrm['rootDia']
    
    cond_VmAm = [
        (RhDFrm['dbh_cm']/100 > RhDFrm['water_depth']),
        (RhDFrm['dbh_cm']/100 <= RhDFrm['water_depth'])
        ]
    try:
        answer_Vm = [
            0,
            np.pi/4*(RhDFrm['dbh_cm']/100)**2*RhDFrm['water_depth'] #m^3 
            ]
            
        RhDFrm['Vm_stem'] = np.select(cond_VmAm, answer_Vm, default='NA').astype(float)
    except:
        RhDFrm['Vm_stem'] = np.nan
    
    try:
        answer_Am = [
            0,
            RhDFrm['dbh_cm']/100*RhDFrm['water_depth'] #m^3 
            ]
            
        RhDFrm['Am_stem'] = np.select(cond_VmAm, answer_Am, default='NA').astype(float)
    except:
        RhDFrm['Am_stem'] = np.nan
    
    
    RhDFrm['Vm_rhiz'] = RhDFrm['Vol_Root'] + RhDFrm['Vm_stem']
    RhDFrm['Am_rhiz'] = RhDFrm['Area_Root'] + RhDFrm['Am_stem']
    
    ## Groupby and sum of the total volume-area
    avicennia_sum = avicennia.groupby(by='cell_number')[['Vm_avc','Am_avc']].sum()
    rhizopora_sum = RhDFrm.groupby(by='cell_number')[['Vm_rhiz','Am_rhiz']].sum()
    
    veg_contribute = pd.concat([avicennia_sum, rhizopora_sum], axis=1)
    veg_contribute['V_total'] = veg_contribute[["Vm_avc", "Vm_rhiz"]].sum(axis=1)
    veg_contribute['A_total'] = veg_contribute[["Am_avc", "Am_rhiz"]].sum(axis=1)
    
    ## Place the Veg Contribute to the same DataFrame with ba and hs
    veg_ba_hs = pd.concat([lookup_hs, lookup_area, veg_contribute], axis=1)
    e_drag = 5 # set to 5m to obtain realistic value for Cd
    veg_ba_hs['V'] = veg_ba_hs['cell_area']*veg_ba_hs['water_depth']
    veg_ba_hs['L'] = (veg_ba_hs['V'] - veg_ba_hs['V_total']) / veg_ba_hs['A_total']
    veg_ba_hs['Cdcalc'] = Cd_no + (e_drag/ veg_ba_hs['L'])
    veg_ba_hs['Cdcalc'] = veg_ba_hs['Cdcalc'].fillna(Cd_no)
    
    ## Place the calculated Cdcalc
    drag_coeff = veg_ba_hs['Cdcalc'].values
    
    # ### Get the calculation per cell number
    # avicennia['V'] = avicennia['cell_area']*avicennia['water_depth']
    # RhDFrm['V'] = RhDFrm['cell_area']*RhDFrm['water_depth']
    
    # # Characteristic Length (m)
    # avicennia['L'] = (avicennia['V'] - avicennia['Vm_avc']) / avicennia['Am_avc']
    # RhDFrm['L'] = (RhDFrm['V'] - RhDFrm['Vm_rhiz']) / RhDFrm['Am_rhiz']
    
    # e_drag = 5 # set to 5m to obtain realistic value for Cd
    
    # # avicennia['Cdcalc'] = Cd_no + (e_drag/ avicennia['L'])
    # avicennia['Cdcalc'] = np.where(avicennia['water_depth'] == 0, np.nan, 
    #                                Cd_no + (e_drag/ avicennia['L']))
    # # RhDFrm['Cdcalc'] = Cd_no + (e_drag/ RhDFrm['L'])
    # RhDFrm['Cdcalc'] = np.where(RhDFrm['water_depth'] == 0, np.nan,
    #                             Cd_no + (e_drag/ RhDFrm['L']))
    
    # ### Groupby and sum
    # avc_grouped = avicennia.groupby(by='cell_number')['Cdcalc'].sum()
    # rzp_grouped = RhDFrm.groupby(by='cell_number')['Cdcalc'].sum()
    
    # avc_grouped_filt = avc_grouped.dropna(axis = 0, how = 'all').replace(0, Cd_no)
    # rzp_grouped_filt = rzp_grouped.dropna(axis = 0, how = 'all').replace(0, Cd_no)
    
    # ## No Veg ## masih belum bener
    
    # no_veg_cd = pd.Series(index_veg_cel, name='Cdcalc')\
    #     .reset_index().rename(columns={'index': 'cell_number'}).set_index('cell_number')
    # no_veg_cd_filt = no_veg_cd.loc[no_veg_cd["Cdcalc"] == 0 ].replace(0,Cd_no)
    # no_veg_cd_filt = no_veg_cd_filt.squeeze()
    
    # ## Append All
    # df_ar = avc_grouped_filt.append(rzp_grouped_filt, ignore_index=False)
    # df_ar_no = df_ar.append(no_veg_cd_filt, ignore_index=False).sort_index()
        
    # drag_coeff = df_ar_no.values
    
    return drag_coeff

def prepCalcDraginVect(model_dfm, read_data):
    read_data_copy = read_data.copy()

    ### Prepare for each species
    ## Avicennia
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen

    avicennia = (read_data_copy.loc[read_data_copy['Species']=='Avicennia_marina']).copy() # for avicennia
    # avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    # Volume of the aerial roots
    avicennia['N_pneu'] = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-avicennia['dbh_cm']))))

    ## Rhizopora
    RhDFrm = (read_data_copy.loc[read_data_copy['Species']=='Rhizopora_apiculata']).copy() # for Rhizopora
    # RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)

    RhDFrm['HRmax'] = ((7.56*RhDFrm['dbh_cm']/100) + 0.5).astype(float)
    RhDFrm['N'] = np.ceil(3.15*RhDFrm['HRmax']**2 + 5.3*RhDFrm['HRmax'] + 0.114).astype(int)
    RhDFrm['HR_interval'] = (RhDFrm['HRmax'] -0.25)/(RhDFrm['N']-1)
    # create linspace of root height as list
    RhDFrm['HR_i'] = RhDFrm.apply(lambda x: np.linspace(0.25,  x['HRmax'], x['N']), axis=1)
    
    return avicennia, RhDFrm

def calcDraginVect_fromPrep(avicennia_pr, RhDFrm_pr, model_dfm):
    
    avicennia = avicennia_pr.copy()
    RhDFrm = RhDFrm_pr.copy()

    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation

    d_pneu = 0.01 # m
    h_pneu = 0.15 # m

    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}

    lookup_area = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'cell_area':ba,
                           })
    lookup_hs = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'water_depth':hs,
                           })
    ###
    val_av = np.vstack(avicennia['cell_number'].values)
    bound_hs = lookup_hs['cell_number'].values == val_av
    avicennia['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])

    val_rz = np.vstack(RhDFrm['cell_number'].values)
    bound_hs = lookup_hs['cell_number'].values == val_rz
    RhDFrm['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])

    # drop based on 0 water level
    avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)

    ## Avicennia
    # if function to adjust the h_pneu for Avicennia
    # this equation is from Du, Qin, et al. 2021
    cond_avc = [
        (avicennia['water_depth'].values < h_pneu),
        (avicennia['water_depth'].values >= h_pneu)
        ]
    val_avc = [
        np.pi*d_pneu*avicennia['water_depth'].values/12,
        np.pi*d_pneu*h_pneu/12
        ]

    avicennia['Vm_pneu'] = np.select(cond_avc, val_avc, default='NA').astype(float)

    avicennia['Vm_pneu_total'] = avicennia['Vm_pneu']*avicennia['N_pneu']
    avicennia['Am_pneu_total'] = d_pneu*h_pneu*avicennia['N_pneu']

    # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
    # height is adjusted from the water depth
    avicennia['Vm_trunk'] = np.pi/4*(avicennia['dbh_cm']/100)**2*avicennia['water_depth']#m^3 
    avicennia['Am_trunk'] = (avicennia['dbh_cm']/100)*avicennia['water_depth'] #m^3 

    # sum of the obstacle volume
    avicennia['Vm_avc'] = avicennia['Vm_pneu_total'] + avicennia['Vm_trunk'] #m^3
    # sum of the obstacle area
    avicennia['Am_avc'] = avicennia['Am_pneu_total'] + avicennia['Am_trunk']#m^2  

    ## Rhizopora
    # explode the DataFrame
    RhDFrm_ex = RhDFrm.explode('HR_i')
    RhDFrm_ex['HR_i'] = RhDFrm_ex['HR_i'].astype(float)
    # calculate as normally
    RhDFrm_ex['x2'] = (1.88*RhDFrm_ex['HR_i'])-0.432
    RhDFrm_ex['rootDia'] = (0.04*RhDFrm_ex['dbh_cm']/100) + (0.005*RhDFrm_ex['HR_i']) + 0.024
    RhDFrm_ex['thetaRoot'] = 20.88*RhDFrm_ex['HR_i'] - 46.87 # in degrees
    RhDFrm_ex['b'] = np.where(RhDFrm_ex['HR_i'] > 2.24,
        0,
        np.tan((20.88*RhDFrm_ex['HR_i'] - 46.87)*np.pi/180)
        )
    RhDFrm_ex['a'] = np.where(RhDFrm_ex['HR_i'] > 2.24, 
             -(RhDFrm_ex['HR_i']/RhDFrm_ex['x2']**2), 
             -((RhDFrm_ex['b']*RhDFrm_ex['x2']+RhDFrm_ex['HR_i'])/RhDFrm_ex['x2']**2)
             )
    RhDFrm_ex['x1'] = np.where(RhDFrm_ex['water_depth'] < RhDFrm_ex['HR_i'],
        (-RhDFrm_ex['b']-np.sqrt(RhDFrm_ex['b']**2-(4*RhDFrm_ex['a']*(RhDFrm_ex['HR_i']-RhDFrm_ex['water_depth']))))/(2*RhDFrm_ex['a']),
        0
        )
    try:
        calc_L_array = RhDFrm_ex[["a", "b", "x1", "x2"]].to_numpy()
        def calc_L(x):
            g = lambda f, a,b: np.sqrt(1+((2*a*f)+b)**2)
            return quad(g, x[2],x[3], args=(x[0],x[1]))[0]
        RhDFrm_ex['L'] = np.apply_along_axis(calc_L, axis=1, arr=calc_L_array)
    
        RhDFrm_ex['VR_i'] = RhDFrm_ex['L']*(RhDFrm_ex['rootDia']/2)**2
        RhDFrm_ex['AR_i'] = RhDFrm_ex['L']*RhDFrm_ex['rootDia']
        RhDFrm_grouped = RhDFrm_ex.groupby(by='Unnamed: 0')[['VR_i','AR_i']].sum().to_numpy()
    
        #add the calculated V and A to main DataFrame
        RhDFrm['Vol_Root'] = RhDFrm_grouped[:,0].tolist()
        RhDFrm['Area_Root'] = RhDFrm_grouped[:,1].tolist()
    except:
        print('initiate run all water level 0m')
        
    
        cond_VmAm = [
            (RhDFrm['dbh_cm']/100 > RhDFrm['water_depth']),
            (RhDFrm['dbh_cm']/100 <= RhDFrm['water_depth'])
            ]
        try:
            answer_Vm = [
                0,
                np.pi/4*(RhDFrm['dbh_cm']/100)**2*RhDFrm['water_depth'] #m^3 
                ]
                
            RhDFrm['Vm_stem'] = np.select(cond_VmAm, answer_Vm, default='NA').astype(float)
        except:
            RhDFrm['Vm_stem'] = np.nan
    
        try:
            answer_Am = [
                0,
                RhDFrm['dbh_cm']/100*RhDFrm['water_depth'] #m^3 
                ]
                
            RhDFrm['Am_stem'] = np.select(cond_VmAm, answer_Am, default='NA').astype(float)
        except:
            RhDFrm['Am_stem'] = np.nan
    
    try:
        RhDFrm['Vm_rhiz'] = RhDFrm['Vol_Root'] + RhDFrm['Vm_stem']
        RhDFrm['Am_rhiz'] = RhDFrm['Area_Root'] + RhDFrm['Am_stem']
    
    except: # this is to compensate all 0m water during initial run
        RhDFrm['Vm_rhiz'] = RhDFrm['Vm_stem']
        RhDFrm['Am_rhiz'] = RhDFrm['Am_stem']    

    ## Groupby and sum of the total volume-area
    avicennia_sum = avicennia.groupby(by='cell_number')[['Vm_avc','Am_avc']].sum()
    rhizopora_sum = RhDFrm.groupby(by='cell_number')[['Vm_rhiz','Am_rhiz']].sum()

    veg_contribute = pd.concat([avicennia_sum, rhizopora_sum], axis=1)
    veg_contribute['V_total'] = veg_contribute[["Vm_avc", "Vm_rhiz"]].sum(axis=1)
    veg_contribute['A_total'] = veg_contribute[["Am_avc", "Am_rhiz"]].sum(axis=1)

    ## Place the Veg Contribute to the same DataFrame with ba and hs
    veg_ba_hs = pd.concat([lookup_hs, lookup_area, veg_contribute], axis=1)
    e_drag = 5 # set to 5m to obtain realistic value for Cd
    veg_ba_hs['V'] = veg_ba_hs['cell_area']*veg_ba_hs['water_depth']
    veg_ba_hs['L'] = (veg_ba_hs['V'] - veg_ba_hs['V_total']) / veg_ba_hs['A_total']
    veg_ba_hs['Cdcalc'] = Cd_no + (e_drag/ veg_ba_hs['L'])
    veg_ba_hs['Cdcalc'] = veg_ba_hs['Cdcalc'].fillna(Cd_no)

    ## Place the calculated Cdcalc
    drag_coeff = veg_ba_hs['Cdcalc'].values
    
    return drag_coeff

def n_prepCalcDraginVect(model_dfm, read_data):
    read_data_copy = read_data.copy()

    ### Prepare for each species
    ## Avicennia
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen

    avicennia = (read_data_copy.loc[read_data_copy['Species']=='Avicennia_marina']).copy() # for avicennia
    # avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    # Volume of the aerial roots
    avicennia['N_pneu'] = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-avicennia['dbh_cm']))))

    ## Rhizopora
    RhDFrm = (read_data_copy.loc[read_data_copy['Species']=='Rhizopora_apiculata']).copy() # for Rhizopora
    # RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)

    RhDFrm['HRmax'] = ((7.56*RhDFrm['dbh_cm']/100) + 0.5).astype(float)
    RhDFrm['N'] = np.ceil(3.15*RhDFrm['HRmax']**2 + 5.3*RhDFrm['HRmax'] + 0.114).astype(int)
    RhDFrm['HR_interval'] = (RhDFrm['HRmax'] -0.25)/(RhDFrm['N']-1)
    
    try:
        # create linspace of root height as list
        RhDFrm['HR_i'] = RhDFrm.apply(lambda x: np.linspace(0.25,  x['HRmax'], x['N']), axis=1)
        
        RhDFrm_ex = RhDFrm.explode('HR_i')
        RhDFrm_ex['HR_i'] = RhDFrm_ex['HR_i'].astype(float)
        # calculate as normally
        RhDFrm_ex['x2'] = (1.88*RhDFrm_ex['HR_i'])-0.432
        RhDFrm_ex['rootDia'] = (0.04*RhDFrm_ex['dbh_cm']/100) + (0.005*RhDFrm_ex['HR_i']) + 0.024
        RhDFrm_ex['thetaRoot'] = 20.88*RhDFrm_ex['HR_i'] - 46.87 # in degrees
        RhDFrm_ex['b'] = np.where(RhDFrm_ex['HR_i'] > 2.24,
            0,
            np.tan((20.88*RhDFrm_ex['HR_i'] - 46.87)*np.pi/180)
            )
        RhDFrm_ex['a'] = np.where(RhDFrm_ex['HR_i'] > 2.24, 
                 -(RhDFrm_ex['HR_i']/RhDFrm_ex['x2']**2), 
                 -((RhDFrm_ex['b']*RhDFrm_ex['x2']+RhDFrm_ex['HR_i'])/RhDFrm_ex['x2']**2)
                 )
    except:
        print('no rhizopora')
        RhDFrm_ex = RhDFrm.copy()
        
    return avicennia, RhDFrm, RhDFrm_ex

def n_n_prepCalcDraginVect(read_data, SYS_APP):
    read_data_copy = read_data.copy()

    ### Prepare for each species
    ## Avicennia
    f_pneu = 0.3 # const in van Maanen
    D05_pneu = 20 # const in van Maanen

    avicennia = (read_data_copy.loc[read_data_copy['Species']=='Avicennia_marina']).copy() # for avicennia
    # avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    # Volume of the aerial roots
    avicennia['N_pneu'] = 10025*(1/(1+np.exp(f_pneu*(D05_pneu-avicennia['dbh_cm']))))

    ## Rhizopora
    RhDFrm = (read_data_copy.loc[read_data_copy['Species']=='Rhizopora_apiculata']).copy() # for Rhizopora
    # RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)
    
    lookup_volume = np.load(os.path.join(SYS_APP,'array_volume_rzp.npy')) 
    lookup_area = np.load(os.path.join(SYS_APP,'array_area_rzp.npy')) 
    dbh_is = np.arange(0.1,100.5,0.5) # dbh is from 0.1cm to 100cm with 0.5cm increment
    wat_depth_is = np.arange(0.01,10.01,0.01) # water depth is from 0.01m to 10.01m with 0.01m increment

    g = interp2d(wat_depth_is,dbh_is,lookup_volume, kind='linear')
    h = interp2d(wat_depth_is,dbh_is,lookup_area, kind='linear')
        
    return avicennia, RhDFrm, g, h

def n_calcDraginVect_fromPrep(avicennia_pr, RhDFrm_pr, RhDFrm_pr_ex, model_dfm, drag_coeff):
    
    avicennia = avicennia_pr.copy()
    RhDFrm = RhDFrm_pr.copy()
    RhDFrm_ex = RhDFrm_pr_ex.copy()

    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation

    d_pneu = 0.01 # m
    h_pneu = 0.15 # m

    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}

    lookup_area = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'cell_area':ba,
                           })
    lookup_hs = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'water_depth':hs,
                           })
    ###
    try:
        val_av = np.vstack(avicennia['cell_number'].values)
        bound_hs = lookup_hs['cell_number'].values == val_av
        avicennia['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
    except: 
        print('no avicennia')
    
    try:
        val_rz = np.vstack(RhDFrm_ex['cell_number'].values)
        bound_hs = lookup_hs['cell_number'].values == val_rz
        RhDFrm_ex['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
        
        val_rzo = np.vstack(RhDFrm['cell_number'].values)
        bound_hs = lookup_hs['cell_number'].values == val_rzo
        RhDFrm['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
    except:
        print('no rhizopora')

    # drop based on 0 water level
    try:
        avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    except:
        print('no avicennia')
    try:
        RhDFrm_ex.drop(RhDFrm_ex.loc[RhDFrm_ex['water_depth']==0].index, inplace=True)
    except:
        print('no rhizopora')

    ## Avicennia
    # if function to adjust the h_pneu for Avicennia
    # this equation is from Du, Qin, et al. 2021
    try:
        cond_avc = [
            (avicennia['water_depth'].values < h_pneu),
            (avicennia['water_depth'].values >= h_pneu)
            ]
        val_avc = [
            np.pi*d_pneu*avicennia['water_depth'].values/12,
            np.pi*d_pneu*h_pneu/12
            ]
    
        avicennia['Vm_pneu'] = np.select(cond_avc, val_avc, default='NA').astype(float)
    
        avicennia['Vm_pneu_total'] = avicennia['Vm_pneu']*avicennia['N_pneu']
        avicennia['Am_pneu_total'] = d_pneu*h_pneu*avicennia['N_pneu']
    
        # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
        # height is adjusted from the water depth
        avicennia['Vm_trunk'] = np.pi/4*(avicennia['dbh_cm']/100)**2*avicennia['water_depth']#m^3 
        avicennia['Am_trunk'] = (avicennia['dbh_cm']/100)*avicennia['water_depth'] #m^3 
    
        # sum of the obstacle volume
        avicennia['Vm_avc'] = avicennia['Vm_pneu_total'] + avicennia['Vm_trunk'] #m^3
        # sum of the obstacle area
        avicennia['Am_avc'] = avicennia['Am_pneu_total'] + avicennia['Am_trunk']#m^2  
    except:
        print('no avicennia')

    ## Rhizopora
    try:
        RhDFrm_ex['x1'] = np.where(RhDFrm_ex['water_depth'] < RhDFrm_ex['HR_i'],
            (-RhDFrm_ex['b']-np.sqrt(RhDFrm_ex['b']**2-(4*RhDFrm_ex['a']*(RhDFrm_ex['HR_i']-RhDFrm_ex['water_depth']))))/(2*RhDFrm_ex['a']),
            0
            )
        try:
            calc_L_array = RhDFrm_ex[["a", "b", "x1", "x2"]].to_numpy()
            def calc_L(x):
                g = lambda f, a,b: np.sqrt(1+((2*a*f)+b)**2)
                return quad(g, x[2],x[3], args=(x[0],x[1]))[0]
            RhDFrm_ex['L'] = np.apply_along_axis(calc_L, axis=1, arr=calc_L_array)
        
            RhDFrm_ex['VR_i'] = RhDFrm_ex['L']*(RhDFrm_ex['rootDia']/2)**2
            RhDFrm_ex['AR_i'] = RhDFrm_ex['L']*RhDFrm_ex['rootDia']
            RhDFrm_grouped = RhDFrm_ex.groupby(by='Unnamed: 0')[['VR_i','AR_i']].sum().to_numpy()
        
            #add the calculated V and A to main DataFrame
            RhDFrm['Vol_Root'] = RhDFrm_grouped[:,0].tolist()
            RhDFrm['Area_Root'] = RhDFrm_grouped[:,1].tolist()
        except:
            print('initiate run all water level 0m')
            
        
            cond_VmAm = [
                (RhDFrm['dbh_cm']/100 > RhDFrm['water_depth']),
                (RhDFrm['dbh_cm']/100 <= RhDFrm['water_depth'])
                ]
            try:
                answer_Vm = [
                    0,
                    np.pi/4*(RhDFrm['dbh_cm']/100)**2*RhDFrm['water_depth'] #m^3 
                    ]
                    
                RhDFrm['Vm_stem'] = np.select(cond_VmAm, answer_Vm, default='NA').astype(float)
            except:
                RhDFrm['Vm_stem'] = np.nan
        
            try:
                answer_Am = [
                    0,
                    RhDFrm['dbh_cm']/100*RhDFrm['water_depth'] #m^3 
                    ]
                    
                RhDFrm['Am_stem'] = np.select(cond_VmAm, answer_Am, default='NA').astype(float)
            except:
                RhDFrm['Am_stem'] = np.nan
        
        try:
            RhDFrm['Vm_rhiz'] = RhDFrm['Vol_Root'] + RhDFrm['Vm_stem']
            RhDFrm['Am_rhiz'] = RhDFrm['Area_Root'] + RhDFrm['Am_stem']
        
        except: # this is to compensate all 0m water during initial run
            RhDFrm['Vm_rhiz'] = RhDFrm['Vm_stem']
            RhDFrm['Am_rhiz'] = RhDFrm['Am_stem']    
    except:
        print('no rhizopora')

    ## Groupby and sum of the total volume-area
    try:
        avicennia_sum = avicennia.groupby(by='cell_number')[['Vm_avc','Am_avc']].sum()
    except:
        print('no avicennia')
    try:
        rhizopora_sum = RhDFrm.groupby(by='cell_number')[['Vm_rhiz','Am_rhiz']].sum()
        veg_contribute = pd.concat([avicennia_sum, rhizopora_sum], axis=1)
    except:
        print('no rhizopora')

    try:
        veg_contribute['V_total'] = veg_contribute[["Vm_avc", "Vm_rhiz"]].sum(axis=1)
        veg_contribute['A_total'] = veg_contribute[["Am_avc", "Am_rhiz"]].sum(axis=1)
    except:
        veg_contribute = pd.DataFrame()
        try:
            veg_contribute['V_total'] = avicennia_sum['Vm_avc']
            veg_contribute['A_total'] = avicennia_sum['Am_avc']
        except:
            veg_contribute['V_total'] = RhDFrm['Vm_rhiz']
            veg_contribute['A_total'] = RhDFrm['Am_rhiz']

    ## Place the Veg Contribute to the same DataFrame with ba and hs
    veg_ba_hs = pd.concat([lookup_hs, lookup_area, veg_contribute], axis=1)
    e_drag = 5 # set to 5m to obtain realistic value for Cd
    veg_ba_hs['V'] = veg_ba_hs['cell_area']*veg_ba_hs['water_depth']
    veg_ba_hs['L'] = np.where(veg_ba_hs['V'] < veg_ba_hs['V_total'],
                    veg_ba_hs['V'] / veg_ba_hs['A_total'],
                    (veg_ba_hs['V'] - veg_ba_hs['V_total']) / veg_ba_hs['A_total'])
    veg_ba_hs['Cdcalc'] = Cd_no + (e_drag/ veg_ba_hs['L'])
    veg_ba_hs['Cdcalc'] = veg_ba_hs['Cdcalc'].fillna(Cd_no)
    drag_prev = pd.DataFrame(drag_coeff, columns = ['drag_prev'])
    veg_ba_hs['Cdcalc'] = np.where(veg_ba_hs['water_depth'] == 0, 
                                   drag_prev['drag_prev'], veg_ba_hs['Cdcalc'])

    ## Place the calculated Cdcalc
    drag_coeff = veg_ba_hs['Cdcalc'].values
    
    return drag_coeff


def n_n_calcDraginVect_fromPrep(avicennia_pr, RhDFrm_pr, vol_func, area_func, model_dfm, drag_coeff):
    
    avicennia = avicennia_pr.copy()
    RhDFrm = RhDFrm_pr.copy()
    # RhDFrm_ex = RhDFrm_pr_ex.copy()

    Cd_no = 0.005 # assume drag coefficient without the presence of vegetation

    d_pneu = 0.01 # m
    h_pneu = 0.15 # m

    ba = model_dfm.get_var('ba') #surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
    hs = model_dfm.get_var('hs') #water depth at the end of timestep {"location": "face", "shape": ["ndx"]}

    lookup_area = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'cell_area':ba,
                           })
    lookup_hs = pd.DataFrame({'cell_number':list(range(len(ba))),
                           'water_depth':hs,
                           })
    ###
    try:
        val_av = np.vstack(avicennia['cell_number'].values)
        bound_hs = lookup_hs['cell_number'].values == val_av
        avicennia['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
    except: 
        print('assign water depth: no avicennia')
    
    try:
        # val_rz = np.vstack(RhDFrm_ex['cell_number'].values)
        # bound_hs = lookup_hs['cell_number'].values == val_rz
        # RhDFrm_ex['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
        
        val_rzo = np.vstack(RhDFrm['cell_number'].values)
        bound_hs = lookup_hs['cell_number'].values == val_rzo
        RhDFrm['water_depth'] = np.dot(bound_hs, lookup_hs['water_depth'])
    except:
        print('assign water depth: no rhizopora')

    # drop based on 0 water level
    try:
        avicennia.drop(avicennia.loc[avicennia['water_depth']==0].index, inplace=True)
    except:
        print('drop water depth 0m: no avicennia')
    try:
        RhDFrm.drop(RhDFrm.loc[RhDFrm['water_depth']==0].index, inplace=True)
    except:
        print('drop water depth 0m:no rhizopora')

    ## Avicennia
    # if function to adjust the h_pneu for Avicennia
    # this equation is from Du, Qin, et al. 2021
    try:
        if avicennia.size == 0:
            print('calc Vm-Am: no avicennia due to water depth 0m')
        else:
            cond_avc = [
                (avicennia['water_depth'].values < h_pneu),
                (avicennia['water_depth'].values >= h_pneu)
                ]
            val_avc = [
                np.pi*d_pneu*avicennia['water_depth'].values/12,
                np.pi*d_pneu*h_pneu/12
                ]
        
            avicennia['Vm_pneu'] = np.select(cond_avc, val_avc, default='NA').astype(float)
        
            avicennia['Vm_pneu_total'] = avicennia['Vm_pneu']*avicennia['N_pneu']
            avicennia['Am_pneu_total'] = d_pneu*h_pneu*avicennia['N_pneu']
        
            # Calculate the d_0 from height of the tree, wait for Uwe Grueter Confirmation
            # height is adjusted from the water depth
            avicennia['Vm_trunk'] = np.pi/4*(avicennia['dbh_cm']/100)**2*avicennia['water_depth']#m^3 
            avicennia['Am_trunk'] = (avicennia['dbh_cm']/100)*avicennia['water_depth'] #m^3 
        
            # sum of the obstacle volume
            avicennia['Vm_avc'] = avicennia['Vm_pneu_total'] + avicennia['Vm_trunk'] #m^3
            # sum of the obstacle area
            avicennia['Am_avc'] = avicennia['Am_pneu_total'] + avicennia['Am_trunk']#m^2  
    except:
        print('calc Vm-Am: no avicennia')

    ## Rhizopora
    try:           
        # RhDFrm['Vm_rhiz'] = g_func(RhDFrm['water_depth'], RhDFrm['dbh_cm']) #water depth(m), dbh(cm)
        
        data_array = RhDFrm[["water_depth", "dbh_cm"]].to_numpy()
        def calc_vol(val):
            return  vol_func(val[0], val[1])
        def calc_are(val):
            return  area_func(val[0], val[1])
        
        RhDFrm['Vm_rhiz'] = np.apply_along_axis(calc_vol, axis=1, arr=data_array)
        RhDFrm['Am_rhiz'] = np.apply_along_axis(calc_are, axis=1, arr=data_array)
        
        # RhDFrm_ex['x1'] = np.where(RhDFrm_ex['water_depth'] < RhDFrm_ex['HR_i'],
        #     (-RhDFrm_ex['b']-np.sqrt(RhDFrm_ex['b']**2-(4*RhDFrm_ex['a']*(RhDFrm_ex['HR_i']-RhDFrm_ex['water_depth']))))/(2*RhDFrm_ex['a']),
        #     0
        #     )
        # try:
        #     calc_L_array = RhDFrm_ex[["a", "b", "x1", "x2"]].to_numpy()
        #     def calc_L(x):
        #         g = lambda f, a,b: np.sqrt(1+((2*a*f)+b)**2)
        #         return quad(g, x[2],x[3], args=(x[0],x[1]))[0]
        #     RhDFrm_ex['L'] = np.apply_along_axis(calc_L, axis=1, arr=calc_L_array)
        
        #     RhDFrm_ex['VR_i'] = RhDFrm_ex['L']*(RhDFrm_ex['rootDia']/2)**2
        #     RhDFrm_ex['AR_i'] = RhDFrm_ex['L']*RhDFrm_ex['rootDia']
        #     RhDFrm_grouped = RhDFrm_ex.groupby(by='Unnamed: 0')[['VR_i','AR_i']].sum().to_numpy()
        
        #     #add the calculated V and A to main DataFrame
        #     RhDFrm['Vol_Root'] = RhDFrm_grouped[:,0].tolist()
        #     RhDFrm['Area_Root'] = RhDFrm_grouped[:,1].tolist()
        # except:
        #     print('initiate run all water level 0m')
            
        
        #     cond_VmAm = [
        #         (RhDFrm['dbh_cm']/100 > RhDFrm['water_depth']),
        #         (RhDFrm['dbh_cm']/100 <= RhDFrm['water_depth'])
        #         ]
        #     try:
        #         answer_Vm = [
        #             0,
        #             np.pi/4*(RhDFrm['dbh_cm']/100)**2*RhDFrm['water_depth'] #m^3 
        #             ]
                    
        #         RhDFrm['Vm_stem'] = np.select(cond_VmAm, answer_Vm, default='NA').astype(float)
        #     except:
        #         RhDFrm['Vm_stem'] = np.nan
        
        #     try:
        #         answer_Am = [
        #             0,
        #             RhDFrm['dbh_cm']/100*RhDFrm['water_depth'] #m^3 
        #             ]
                    
        #         RhDFrm['Am_stem'] = np.select(cond_VmAm, answer_Am, default='NA').astype(float)
        #     except:
        #         RhDFrm['Am_stem'] = np.nan
        
        # try:
        #     RhDFrm['Vm_rhiz'] = RhDFrm['Vol_Root'] + RhDFrm['Vm_stem']
        #     RhDFrm['Am_rhiz'] = RhDFrm['Area_Root'] + RhDFrm['Am_stem']
        
        # except: # this is to compensate all 0m water during initial run
        #     RhDFrm['Vm_rhiz'] = RhDFrm['Vm_stem']
        #     RhDFrm['Am_rhiz'] = RhDFrm['Am_stem']    
    except:
        print('calc Vm-Am: no rhizopora')

    ## Groupby and sum of the total volume-area
    try:
        avicennia_sum = avicennia.groupby(by='cell_number')[['Vm_avc','Am_avc']].sum()
    except:
        print('calc groupby: no avicennia')
    try:
        rhizopora_sum = RhDFrm.groupby(by='cell_number')[['Vm_rhiz','Am_rhiz']].sum()
        veg_contribute = pd.concat([avicennia_sum, rhizopora_sum], axis=1)
    except:
        print('calc groupby: no rhizopora')
        
    try:
        try:
            veg_contribute['V_total'] = veg_contribute[["Vm_avc", "Vm_rhiz"]].sum(axis=1)
            veg_contribute['A_total'] = veg_contribute[["Am_avc", "Am_rhiz"]].sum(axis=1)
        except:
            veg_contribute = pd.DataFrame()
            try:
                veg_contribute['V_total'] = avicennia_sum['Vm_avc']
                veg_contribute['A_total'] = avicennia_sum['Am_avc']
            except:
                veg_contribute['V_total'] = RhDFrm['Vm_rhiz']
                veg_contribute['A_total'] = RhDFrm['Am_rhiz']
    except:
        print('no rhizopora and no avicennia')
        total_va = np.array([np.zeros(len(lookup_hs)), np.zeros(len(lookup_hs))]).T
        total_va[:] = np.nan
        veg_contribute = pd.DataFrame(data = total_va ,columns=['V_total', 'A_total'])

    ## Place the Veg Contribute to the same DataFrame with ba and hs
    veg_ba_hs = pd.concat([lookup_hs, lookup_area, veg_contribute], axis=1)
    e_drag = 5 # set to 5m to obtain realistic value for Cd
    veg_ba_hs['V'] = veg_ba_hs['cell_area']*veg_ba_hs['water_depth']
    veg_ba_hs['L'] = np.where(veg_ba_hs['V'] < veg_ba_hs['V_total'],
                    veg_ba_hs['V'] / veg_ba_hs['A_total'],
                    (veg_ba_hs['V'] - veg_ba_hs['V_total']) / veg_ba_hs['A_total'])
    veg_ba_hs['Cdcalc'] = Cd_no + (e_drag/ veg_ba_hs['L'])
    veg_ba_hs['Cdcalc'] = veg_ba_hs['Cdcalc'].fillna(Cd_no)
    drag_prev = pd.DataFrame(drag_coeff, columns = ['drag_prev'])
    veg_ba_hs['Cdcalc'] = np.where(veg_ba_hs['water_depth'] == 0, 
                                   drag_prev['drag_prev'], veg_ba_hs['Cdcalc'])
    ## control the Cdcalc in a very shallow environment, use Cd range as in Mazda, 1997 (0.4-10)
    veg_ba_hs['Cdcalc'] = np.where(veg_ba_hs['Cdcalc'] > 10,10, veg_ba_hs['Cdcalc'] )

    ## Place the calculated Cdcalc
    drag_coeff = veg_ba_hs['Cdcalc'].values
    
    return drag_coeff

def calc_woo_avc_rzp(h_wl, blv_wl, model_dfm):
    test_logic = h_wl > blv_wl # previous woo
    # test_logic = h_wl <= blv_wl # true means dry
    def calc_woo_per_cell(test_logic):
        test_pd = pd.DataFrame(test_logic, columns=['data_0'])
        test_pd['avc_1'] = test_pd['data_0'].shift(periods=-1)
        test_pd['avc_2'] = test_pd['data_0'].shift(periods=-2)
        test_pd['sum_avc'] = test_pd['data_0'] + test_pd['avc_1'] + test_pd['avc_2']
        
        test_pd['rzp_3'] = test_pd['data_0'].shift(periods=-3)
        test_pd['rzp_4'] = test_pd['data_0'].shift(periods=-4)
        test_pd['sum_rzp'] = test_pd['sum_avc'] + test_pd['rzp_3'] + test_pd['rzp_4']
        
        test_pd['sol_av'] = np.where(test_pd['sum_avc']==0,1,0) # sum_avc 0 means dry within threshold, therefore score 1
        test_pd['sol_rzp'] = np.where(test_pd['sum_rzp']==0,1,0)
        # test_pd['sol_av'] = np.where(test_pd['sum_avc']==3,1,0) # sum_avc 3 means dry within threshold, therefore score 1
        # test_pd['sol_rzp'] = np.where(test_pd['sum_rzp']==5,1,0)
        
        val_is = np.empty(2)
        val_is[0] = (len(test_pd.index)-test_pd['sol_av'].sum()) / len(test_pd.index)
        val_is[1] = (len(test_pd.index)-test_pd['sol_rzp'].sum()) / len(test_pd.index)
        # val_is[0] = test_pd['sol_av'].sum() / len(test_pd.index)
        # val_is[1] = test_pd['sol_rzp'].sum() / len(test_pd.index)
        
        # return val_avc, val_rzp
        return val_is
    
    val_cell = np.empty([test_logic.shape[0],2])
    bed_is_logic = model_dfm.get_var('bl') >= 0
    for ii in range(val_cell.shape[0]):
        if bed_is_logic[ii] == False:
            val_cell[ii,:] = 0
        else:
            val_cell[ii,:] = calc_woo_per_cell(test_logic[ii])
    
    # bed_is_logic = np.tile((model_dfm.get_var('bl') >= 0), (2,1)).T
    # val_cell_last = np.where(bed_is_logic == True, val_cell, 0)
    
        
    return val_cell

def calc_woo_avc_rzp_hs(hs_wl, model_dfm):
    test_logic = hs_wl <= 0.01 # where water less or equal to 1cm it is dry
    # based on Sousa, et.al(2003) paper regarding propagule dimension
    # test_logic = hs_wl == 0 # where completely dry
    # test_logic = h_wl <= blv_wl # true means dry
    def calc_woo_per_cell(test_logic):
        test_pd = pd.DataFrame(test_logic, columns=['data_0'])
        test_pd['avc_1'] = test_pd['data_0'].shift(periods=-1)
        test_pd['avc_2'] = test_pd['data_0'].shift(periods=-2)
        test_pd['sum_avc'] = test_pd['data_0'] + test_pd['avc_1'] + test_pd['avc_2']
        
        test_pd['rzp_3'] = test_pd['data_0'].shift(periods=-3)
        test_pd['rzp_4'] = test_pd['data_0'].shift(periods=-4)
        test_pd['sum_rzp'] = test_pd['sum_avc'] + test_pd['rzp_3'] + test_pd['rzp_4']
        
        # test_pd['sol_av'] = np.where(test_pd['sum_avc']==0,1,0) # sum_avc 0 means dry within threshold, therefore score 1
        # test_pd['sol_rzp'] = np.where(test_pd['sum_rzp']==0,1,0)
        test_pd['sol_av'] = np.where(test_pd['sum_avc']==3,1,0) # sum_avc 3 means dry within threshold, therefore score 1
        test_pd['sol_rzp'] = np.where(test_pd['sum_rzp']==5,1,0)
        
        val_is = np.empty(2)
        # val_is[0] = (len(test_pd.index)-test_pd['sol_av'].sum()) / len(test_pd.index)
        # val_is[1] = (len(test_pd.index)-test_pd['sol_rzp'].sum()) / len(test_pd.index)
        val_is[0] = test_pd['sol_av'].sum() / len(test_pd.index) # probabillity of survival
        val_is[1] = test_pd['sol_rzp'].sum() / len(test_pd.index)
        
        # return val_avc, val_rzp
        return val_is
    
    val_cell = np.empty([test_logic.shape[0],2])
    bed_is_logic = model_dfm.get_var('bl') >= 0
    for ii in range(val_cell.shape[0]):
        if bed_is_logic[ii] == False:
            val_cell[ii,:] = 0
        else:
            val_cell[ii,:] = calc_woo_per_cell(test_logic[ii])   
        
    return val_cell