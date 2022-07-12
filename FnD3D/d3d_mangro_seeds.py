# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:26:40 2022

@author: sbe002
"""

import numpy as np
import numpy.ma as ma
import pandas as pd
import random
import math
import datetime
from numpy.random import default_rng

def index_veg(model_dfm, xyzw_cell_number, xyzw_nodes, xk, yk, read_data):
    index_veg_cel = np.empty((model_dfm.get_var('Cdvegsp').shape[0],0))
    for row in range(len(xyzw_cell_number)):
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
        if read_data_subset.shape[0] == 0: #if 0 then skip
            index_is = 0
        else:
            index_is = 1
        index_veg_cel = np.append(index_veg_cel, index_is)
    
    return index_veg_cel

def index_veg_cdveg(xzyz_cell_number, ugrid_all, read_data):
    index_veg_cel = np.empty((xzyz_cell_number.shape[0],0))
    data_verts = ugrid_all.verts
    for row in range(len(xzyz_cell_number)):
        position = xzyz_cell_number[row,2].astype(int)
        aa = data_verts[position,:,:]
        # subsetting pandas 
        read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                     (read_data['GeoRefPosX'] <= max(aa[:,0]))]
        read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                            (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
        if read_data_subset.shape[0] == 0: #if 0 then skip
            index_is = 0
        else:
            index_is = 1
        index_veg_cel = np.append(index_veg_cel, index_is)
    
    return index_veg_cel

class seedling_establishment():
    def __init__(self,data_mangrove):
        self.data_mangrove = data_mangrove
        
    def reduction_nutrient_avicenniamarina(self): ## tidak digunakan saat ini
        c1 = -0.5
        c2 = 2.88
        c3 = -1.66
        RNA = 369.46/1000 ## bukan ini nilainya
        f_nut_red = c1 + c2*RNA + c3*RNA**2 ## harusnya saat RNA max = 1
        
        return f_nut_red
        
    def establishment_avicennia(self, salinity_value):
        d = -0.18
        ui = 72
        D = 0.03 # in R. apiculata D=0.006 (Grueters, et.al., 2014), another set: 0.03
        # from Schile, Lisa M., et al. 2017, seed reproduction is 1.076 per 1m^2 crown area
        f_sal_red = 1/(1+np.exp(d*(ui-salinity_value)))
        # self.data_mangrove['Species'] = self.data_mangrove['Species'].astype('str')
        # it is better to directly calculate seedlings productions based on crown surface
        # not limited to age, as it is difficult to measure age on field
        # mostly via dbh as a proxy
        # avicennia_only = self.data_mangrove[(self.data_mangrove['Species'].str.contains("Avicennia_marina")) & 
        #                                     (self.data_mangrove['Age'] >= 4)]
        avicennia_only = self.data_mangrove[(self.data_mangrove['Species'].str.contains("Avicennia_marina"))] 
        
        N = round(f_sal_red*D*avicennia_only['CrownSurfaceArea_m2'])
        N.rename('N seedlings', inplace=True)
        N_seed = pd.merge(N,avicennia_only, left_index=True, right_index=True)
        N_seed.drop(['dbh_cm', 'Height_cm','Age'], inplace=True, axis=1)
        return N_seed
    
# =============================================================================
#     def seedlings_drift(self,data_n, residual,ts): #assume residual current is as pd DataFrame
#         #creatinga a uniformly distributed seedling position
#         # source: https://stackoverflow.com/questions/30564015/how-to-generate-random-points-in-a-circular-distribution
#         xx = []
#         yy = []
#         finalPosX =[]
#         finalPosY = []
#         data_n = data_n.reset_index(drop=True)
#         for n_point in range(data_n.shape[0]):
#             ran_radius = data_n['CrownSurfaceArea_m2'][n_point]
#             pos_x = data_n['GeoRefPosX'][n_point]
#             pos_y = data_n['GeoRefPosY'][n_point]
#             # r, theta = [math.sqrt(random.randint(0,ran_radius))*math.sqrt(ran_radius), 2*math.pi*random.random()]
#             r, theta = [math.sqrt(random.uniform(0,ran_radius))*math.sqrt(ran_radius), 2*math.pi*random.random()]
#             xxi = pos_x + r * math.cos(theta) 
#             yyi = pos_y + r * math.sin(theta)
# 
#             xx = np.append(xx,xxi)
#             yy = np.append(yy,yyi)
#             
#             finalPosX = xx+residual[0]*ts
#             finalPosY = yy+residual[1]*ts
#             # seedsAge = np.zeros(finalPosX.size)
#             
#             # final_position = np.column_stack((finalPosX,finalPosY))
#             # final_position = pd.DataFrame({'seedsPosX': finalPosX, 'seedsPosY': finalPosY, 'seedsAge':seedsAge})
#         final_position = pd.DataFrame({'seedsPosX': finalPosX, 'seedsPosY': finalPosY, 
#                                            'ParentPosX':pos_x, 'ParentPosY':pos_y,
#                                            'CrownSurfaceArea_m2':ran_radius})
#         final_position['Species']=data_n['Species']    
#         
#         return final_position
# =============================================================================



def nn_drift_seed(data_n, residual,ts, n_point, row): #assume residual current is as pd DataFrame
    #creatinga a uniformly distributed seedling position
    # source: https://stackoverflow.com/questions/30564015/how-to-generate-random-points-in-a-circular-distribution
    xx = []
    yy = []
    finalPosX =[]
    finalPosY = []
    data_n = data_n.reset_index(drop=True)
    
    ran_radius = data_n['CrownSurfaceArea_m2'][n_point]
    pos_x = data_n['GeoRefPosX'][n_point]
    pos_y = data_n['GeoRefPosY'][n_point]
    r, theta = [math.sqrt(random.uniform(0,ran_radius))*math.sqrt(ran_radius), 
                2*math.pi*random.random()]
    xxi = pos_x + r * math.cos(theta) 
    yyi = pos_y + r * math.sin(theta)

    xx = np.append(xx,xxi)
    yy = np.append(yy,yyi)
    
    finalPosX = xx+residual[0]*ts
    finalPosY = yy+residual[1]*ts

    final_position = pd.DataFrame({'seedsPosX': finalPosX, 'seedsPosY': finalPosY, 
                                   'ParentPosX':pos_x, 'ParentPosY':pos_y,
                                   'CrownSurfaceArea_m2':ran_radius})
    final_position['Species']=data_n['Species']
    final_position['row']=row
    
    return final_position

def subsetting_cell(xyzw_cell_number, row, xyzw_nodes, xk, yk, read_data):
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
    
    return read_data_subset

def subsetting_cellCdveg(ugrid_all, xzyz_cell_number, row, read_data):
    data_verts = ugrid_all.verts
    position = xzyz_cell_number[row,2].astype(int)
     
    aa = data_verts[position,:,:]
    # subsetting pandas 
    read_data_subset = read_data[(read_data['GeoRefPosX'] >= min(aa[:,0])) & 
                                 (read_data['GeoRefPosX'] <= max(aa[:,0]))]
    read_data_subset = read_data_subset[(read_data_subset['GeoRefPosY'] >= min(aa[:,1])) & 
                                        (read_data_subset['GeoRefPosY'] <= max(aa[:,1]))]
    
    return read_data_subset

# =============================================================================
# def seedling_dispersal(xyzw_cell_number, index_veg_cel, xyzw_nodes, xk, yk, 
#                        read_data, med_sal, residual_is, model_dfm, reftime, cursec):
#     list_of_seeds = []
#     for row in range(len(xyzw_cell_number)):
#         if index_veg_cel[row] == 1:
#             read_data_subset = subsetting_cell(xyzw_cell_number, row, 
#                                                xyzw_nodes, xk, yk, read_data)  
#             seedss = seedling_establishment(read_data_subset)
#             Nn = seedss.establishment_avicennia(med_sal[row])
#             if Nn.size > 0:
#                 seedling_pos = seedss.seedlings_drift(Nn, residual_is[row],
#                                                       model_dfm.get_time_step())
#                 list_of_seeds.append(seedling_pos)
#     
#     seedling_finalpos = pd.concat(list_of_seeds, axis=0)
#     seedling_finalpos['Age'] = datetime.timedelta(days = 0)
#     # modify the column name
#     seedling_finalpos = seedling_finalpos.rename(columns={'seedsPosX':'GeoRefPosX',
#                                                   'seedsPosY':'GeoRefPosY'})
#     seedling_finalpos['Created Time Stamp'] = reftime + cursec
#     seedling_finalpos = seedling_finalpos.reset_index(drop=True)
#     
#     return seedling_finalpos
# =============================================================================

# def seedling_dispersal(xyzw_cell_number, index_veg_cel, xyzw_nodes, xk, yk, 
#                        read_data, med_sal, residual_is, model_dfm, reftime, cursec):
def seedling_dispersal(xzyz_cell_number, index_veg_cel, ugrid_all, read_data, 
                       med_sal, residual_is, model_dfm, reftime, cursec):
    list_of_seeds = [] 
    for row in range(len(xzyz_cell_number)):
    # for row in range(len(xyzw_cell_number)):
        # print(row) # row yang ada veg adalah 78
        if index_veg_cel[row] == 1:
            # print (index_veg_cel[row] == 1, row)
            # read_data_subset = subsetting_cell(xyzw_cell_number, row, 
            #                                    xyzw_nodes, xk, yk, read_data)  
            read_data_subset = subsetting_cellCdveg(ugrid_all, xzyz_cell_number, 
                                                    row, read_data)
            seedss = seedling_establishment(read_data_subset)
            Nn = seedss.establishment_avicennia(med_sal[row])
            Nn = Nn.reset_index(drop=True)
            seedling_pos_lists = []
            for n_point in range(Nn.shape[0]):
                # print(n_point)
                if Nn['N seedlings'][n_point] > 0:
                    for n_seeds in range(Nn['N seedlings'][n_point].astype(int)):
                        # print(row)
                        seedling_at_loop = nn_drift_seed(Nn, residual_is[row],
                                            model_dfm.get_time_step(), n_point,
                                            row)
                        seedling_pos_lists.append(seedling_at_loop)
                    # seedling_final_pos = pd.concat(seedling_pos_lists) # sebelumnya pakai ini tanpa try except
                    try:
                        seedling_final_pos = pd.concat(seedling_pos_lists)
                    except ValueError: 
                        seedling_final_pos = pd.DataFrame()
                        # seedling_final_pos = []
                else:    
                    pass
                    
            try:
                list_of_seeds.append(seedling_final_pos)
            except UnboundLocalError:
                pass

        else:
            pass

    
    # seedling_finalpos = pd.concat(list_of_seeds, axis=0) # sebelumnya pakai ini semua tanpa try except
    # seedling_finalpos['Age'] = datetime.timedelta(days = 0)
    # # modify the column name
    # seedling_finalpos = seedling_finalpos.rename(columns={'seedsPosX':'GeoRefPosX',
    #                                               'seedsPosY':'GeoRefPosY'})
    # seedling_finalpos['Created Time Stamp'] = reftime + cursec
    # seedling_finalpos = seedling_finalpos.reset_index(drop=True)
    try:
        seedling_finalpos = pd.concat(list_of_seeds, axis=0)
        seedling_finalpos['Age'] = datetime.timedelta(days = 0)
        # modify the column name
        seedling_finalpos = seedling_finalpos.rename(columns={'seedsPosX':'GeoRefPosX',
                                                      'seedsPosY':'GeoRefPosY'})
        seedling_finalpos['Created Time Stamp'] = reftime + cursec
        seedling_finalpos = seedling_finalpos.reset_index(drop=True)
    except ValueError:
        seedling_finalpos = pd.DataFrame()
        # seedling_finalpos = []
    
    return seedling_finalpos

# function to calculate residual current   
def calculate_residual(the_res, the_ts, xz):
    x_test = np.empty((len(xz),0))
    for column in range(the_res.shape[1]-1):
        x_sum = the_res[:,column:column+2].sum(axis=1)/the_ts
        x_test = np.append(x_test, np.reshape(x_sum,(len(x_sum),1)), axis=1)
    residual_of = np.median(x_test, axis=1)
    
    return residual_of
def collect_res( model_dfm, res_curr_x, res_curr_y):
    velocity_in_x = model_dfm.get_var('ucx') # velocity x
    velocity_in_y = model_dfm.get_var('ucy') # velocity y
    res_curr_x = np.append(res_curr_x, np.reshape(velocity_in_x,
                            (len(velocity_in_x),1)), axis=1)
    res_curr_y = np.append(res_curr_y, np.reshape(velocity_in_y,
                            (len(velocity_in_y),1)), axis=1)
    return res_curr_x, res_curr_y

# function to filter the probability of the dispersed seedlings due to WoO

# def seedling_prob(sdlg_fnl_pos, xyzw_cell_number,xyzw_nodes, xk, yk, surv_val):
def seedling_prob(sdlg_fnl_pos, xzyz_cell_number, ugrid_all, surv_val):
    list_of_seeds = [] 
    # for row in range(len(xyzw_cell_number)):
    for row in range(len(xzyz_cell_number)):
        # print(row) # row yang ada veg adalah 78
        # div_of_seeds = []
        # if index_veg_cel[row] == 1:
            # print (index_veg_cel[row] == 1, row)
        # read_data_subset = subsetting_cell(xyzw_cell_number, row, 
        #                                    xyzw_nodes, xk, yk, sdlg_fnl_pos)
        read_data_subset = subsetting_cellCdveg(ugrid_all, xzyz_cell_number, row, sdlg_fnl_pos)
        surv_val_in_row = surv_val[row]
        
        if surv_val_in_row <= 0.3:
            num_of_seeds = read_data_subset.copy()
            div_num_seeds = round(num_of_seeds.shape[0]/1.5)
            if div_num_seeds == 0:
                num_of_seeds = pd.DataFrame()
            else:
                size_is = num_of_seeds.shape[0]-div_num_seeds
                arr_indices_top_drop = default_rng().choice(num_of_seeds.index, size=size_is, replace=False)
                after_calc = num_of_seeds.drop(index=arr_indices_top_drop)
                num_of_seeds = after_calc
        elif surv_val_in_row <= 0.5:
            num_of_seeds = read_data_subset.copy()
            div_num_seeds = round(num_of_seeds.shape[0]/2)
            if div_num_seeds == 0:
                num_of_seeds = pd.DataFrame()
            else:
                size_is = num_of_seeds.shape[0]-div_num_seeds
                arr_indices_top_drop = default_rng().choice(num_of_seeds.index, size=size_is, replace=False)
                after_calc = num_of_seeds.drop(index=arr_indices_top_drop)
                num_of_seeds = after_calc
        elif surv_val_in_row <= 0.7:
            num_of_seeds = read_data_subset.copy()
            div_num_seeds = round(num_of_seeds.shape[0]/3)
            if div_num_seeds == 0:
                num_of_seeds = pd.DataFrame()
            else:
                size_is = num_of_seeds.shape[0]-div_num_seeds
                arr_indices_top_drop = default_rng().choice(num_of_seeds.index, size=size_is, replace=False)
                after_calc = num_of_seeds.drop(index=arr_indices_top_drop)
                num_of_seeds = after_calc
        else:
            num_of_seeds = read_data_subset.copy()
    
        try:
            
            list_of_seeds.append(num_of_seeds)
        except:
            pass
    
    seedling_after_woo = pd.concat(list_of_seeds, axis=0)
    seedling_after_woo = seedling_after_woo.reset_index(drop=True)
    return seedling_after_woo
            
def elim_seeds_surv(sdlg_fnl_pos, xzyz_cell_number, ugrid_all, surv_val):
    list_of_seeds = []
    
    for row in range(len(xzyz_cell_number)):
        surv_val_in_row = surv_val[row]
        if surv_val_in_row > 0:
            read_data_subset = subsetting_cellCdveg(ugrid_all, xzyz_cell_number, row, sdlg_fnl_pos)
            # seeds_not_zero = read_data_subset.copy()
            list_of_seeds.append(read_data_subset)
    
    seedling_after_filter = pd.concat(list_of_seeds, axis=0)
    seedling_after_filter = seedling_after_filter.reset_index(drop=True)
    return seedling_after_filter

def range_seed(bed_level, limit_seed, LLWL, no_data_val, xz, yz):
    # get the last column of bed_level
    bed_level_last = bed_level[:,-1]
    bed_level = np.where(bed_level_last < LLWL, no_data_val, bed_level_last) # less than LLWL will be -999
    # arrange the xyz
    bed_level_val_raster = np.column_stack((xz,yz,bed_level))
    ## check and filter the bed level if it is within range that you want
    # filter x value between range x from bed_level_raster
    bed_level_check = bed_level_val_raster[ (limit_seed[0] <= bed_level_val_raster[:,0]) 
                                           & (bed_level_val_raster[:,0] <= limit_seed[1]) ] 
    # filter y value between range y from bed_level_check
    bed_level_check = bed_level_check[ (limit_seed[2] <= bed_level_check[:,1]) 
                                      & (bed_level_check[:,1] <= limit_seed[3])]                       
    # filter array for no_data_val
    bed_level_check = bed_level_check[bed_level_check[:,2] != no_data_val]                                                                      
    
    return bed_level_check