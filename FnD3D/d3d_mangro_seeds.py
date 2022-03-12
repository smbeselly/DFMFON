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

class seedling_establishment():
    def __init__(self,data_mangrove):
        self.data_mangrove = data_mangrove
        
    def reduction_nutrient_avicenniamarina(self):
        c1 = -0.5
        c2 = 2.88
        c3 = -1.66
        RNA = 369.46/1000 ## bukan ini nilainya
        f_nut_red = c1 + c2*RNA + c3*RNA**2 ## harusnya saat RNA max = 1
        
        return f_nut_red
        
    def establishment_avicennia(self, salinity_value):
        d = -0.18
        ui = 72
        D = 0.03 # in R. apiculata D=0.006 (Grueters, et.al., 2014)
        f_sal_red = 1/(1+np.exp(d*(ui-salinity_value)))
        # self.data_mangrove['Species'] = self.data_mangrove['Species'].astype('str')
        avicennia_only = self.data_mangrove[(self.data_mangrove['Species'].str.contains("Avicennia_marina")) & 
                                            (self.data_mangrove['Age'] >= 4)]        
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

def seedling_dispersal(xyzw_cell_number, index_veg_cel, xyzw_nodes, xk, yk, 
                       read_data, med_sal, residual_is, model_dfm, reftime, cursec):
    list_of_seeds = [] 
    for row in range(len(xyzw_cell_number)):
        # print(row) # row yang ada veg adalah 78
        if index_veg_cel[row] == 1:
            # print (index_veg_cel[row] == 1, row)
            read_data_subset = subsetting_cell(xyzw_cell_number, row, 
                                               xyzw_nodes, xk, yk, read_data)  
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
                        seedling_final_pos = []
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
        seedling_finalpos = []
    
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
