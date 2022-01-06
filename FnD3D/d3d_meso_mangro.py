# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 19:38:02 2021

@author: sbe002
"""

#%% Import Libraries
import numpy as np
import numpy.ma as ma
import pandas as pd


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

    Cd_calc = Cd_no + (e_drag/L)
    Cd_calc_weighted = cprs_avicennia/ cprs_avicennia.sum() * Cd_calc
    Cd_calc_weighted = Cd_calc_weighted.sum()
    
    return Cd_calc_weighted
