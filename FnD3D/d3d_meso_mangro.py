# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 19:38:02 2021

@author: sbe002
"""

#%% Import Libraries
import numpy as np
import numpy.ma as ma

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