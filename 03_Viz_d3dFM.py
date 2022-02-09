# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 13:41:39 2022

@author: sbe002
"""

from matplotlib import pyplot as plt
import numpy as np
# import geopandas as gpd
import pandas as pd
import numpy as np
import os

PROJ_HOME = os.path.join(r'D:\Git\d3d_meso')
MFON_OUT = os.path.join(PROJ_HOME,'Model-Out','MesoFON')


# Import the last coupling dataset
coup12 = os.path.join(MFON_OUT, 'Compile', 'Coupling_11.txt')
trees_coup12 = pd.read_csv(coup12)

# relative size of the canopy
rr = trees_coup12['Height_cm'].max()/trees_coup12['Height_cm'].min()

# sphere params
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

def canopySphere(trees_coup12, rpar, index, u, v):
    xa = rpar * np.outer(np.cos(u), np.sin(v)) + trees_coup12['GeoRefPosX'][index]
    ya = rpar * np.outer(np.sin(u), np.sin(v)) + trees_coup12['GeoRefPosY'][index]
    za = rpar * np.outer(np.ones(np.size(u)), np.cos(v)) + trees_coup12['Height_cm'][index]
    return xa, ya, za
    
def stemSphere(trees_coup12, index):
    xs = [trees_coup12['GeoRefPosX'][index],trees_coup12['GeoRefPosX'][index]]
    ys = [trees_coup12['GeoRefPosY'][index],trees_coup12['GeoRefPosY'][index]]
    zs = [0,trees_coup12['Height_cm'][index]]
    return xs, ys, zs

def plottingSphere (data, u, v, rr):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for tree in range(data.shape[0]):
        # rpar = data['Height_cm'][tree]/data['Height_cm'].min()
        rpar = data['Height_cm'][tree]/rr
        xa, ya, za = canopySphere(data, rpar, tree, u, v)
        ax.plot_surface(xa, ya, za, color='green')
        
        r1 = data['Height_cm'][tree]/data['Height_cm'].min()
        xs, ys, zs = stemSphere(trees_coup12, tree)
        ax.plot(xs, ys, zs, color='brown',linewidth=r1*2)
        
        ax.view_init(120, 30) # adjust the preview
        ax.set_box_aspect([1,1,1]) 
        print('processing fig ', tree)
        
plottingSphere(trees_coup12, u, v, rr)        

# Tes create the tree as sphere

#plt.rcParams["figure.figsize"] = [7.50, 3.50]
#plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

for tree in range(trees_coup12.shape[0]):
    print(tree)
    rpar = trees_coup12['Height_cm'][tree]/trees_coup12['Height_cm'].min()
    xa, ya, za = canopySphere(trees_coup12, rpar, tree, u, v)
    ax.plot_surface(xa, ya, za, color='green')

for stem in range(trees_coup12.shape[0]):
    print('processing trees ', stem)
    r1 = trees_coup12['Height_cm'][tree]/trees_coup12['Height_cm'].min()
    xs, ys, zs = stemSphere(trees_coup12, stem)
    ax.plot(xs, ys, zs, color='brown',linewidth=r1*2)

ax.set_box_aspect([1,1,1])     
    





# data of the tree position

x, y, z = [1, 1.5], [1, 2.4], [3.4, 1.4]
x,y,z = np.array((x,y,z))

# data for the stem

x2, y2, z2 = [1,1], [1,1], [0, 3.4]
x3, y3, z3 = [1.5,1.5], [2.4,2.4], [0,1.4]

# relative size of the canopy

rr = z[0]/z[1]
r1 = z[0]/rr
r2 = 0.5*r1

# create sphere for the canopy

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
xa = r1 * np.outer(np.cos(u), np.sin(v)) + x[0]
ya = r1 * np.outer(np.sin(u), np.sin(v)) + y[0]
za = r1 * np.outer(np.ones(np.size(u)), np.cos(v)) + z[0]

# 2nd tree

xb = r2 * np.outer(np.cos(u), np.sin(v)) + x[1]
yb = r2 * np.outer(np.sin(u), np.sin(v)) + y[1]
zb = r2 * np.outer(np.ones(np.size(u)), np.cos(v)) + z[1]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.plot_surface(xa, ya, za, color='green')
ax.plot_surface(xb, yb, zb, color='green')

#ax.scatter(x, y, z, c='red', s=z*100)
ax.plot(x2, y2, z2, color='brown',linewidth=r1*2)
ax.plot(x3, y3, z3, color='brown',linewidth=r2*2)

ax.set_box_aspect([1,1,1]) # aspect ratio 1:1:1 in view space
#ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z))) # aspect ratio is 1:1:1 in data space
plt.show()
