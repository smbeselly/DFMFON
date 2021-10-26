# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 13:48:50 2021

@author: sbe002

This is series of slicing examples
"""

#%% Tests
import numpy as np
import numpy.ma as ma

points = np.array([[10,  9], [ 9, 18], [16, 13], [11, 15], [12, 14], [18, 12],
                   [ 2, 14], [ 6, 18], [ 9,  9], [10,  8], [ 6, 17], [ 5,  3],
                   [13, 19], [ 3, 18], [ 8, 17], [ 9,  7], [ 3,  0], [13, 18],
                   [15,  4], [13, 16]])

points2 = np.array([[1. ,2. ,3.], [4. ,5. ,6.]])

#%% Testing
# example = https://note.nkmk.me/en/python-numpy-condition/
aa = np.array([[10,2,np.nan],[41,56,8],[23,5,98],[78,25,1],[235,45,59],
               [65,98,12],[21,89,1],[23,54,7],[41,65,32],[15,62,np.nan]])
#slicing
aax = aa[:,0]
aay = aa[:,1]
aaz = aa[:,2]

print(aa[~np.isnan(aa).any(axis=1)])
print(np.delete(aa, np.where(aa == 12)[0], axis=0))

aa = (aa[~np.isnan(aa).any(axis=1)])
aa= (np.delete(aa, np.where(aa == 12)[0], axis=0))


#%% Example from # source: https://note.nkmk.me/en/python-numpy-nan-remove/
a = np.array([[15,16,78,25],[11, 12, np.nan, 14],[21,np.nan,np.nan,24],[31,32,33,34]])
# Remove all missing values (NaN)
print(np.isnan(a))
print(~np.isnan(a))
print(a[~np.isnan(a)])
print(np.isnan(a).any(axis=1))
print(~np.isnan(a).any(axis=1))
#delete rows with nan
print(a[~np.isnan(a).any(axis=1), :])
print(a[~np.isnan(a).any(axis=1)])
#delete columns with nan
print(~np.isnan(a).any(axis=0))
print(a[:, ~np.isnan(a).any(axis=0)])
# try to use where
print(np.where(np.isnan(a)))
print(a[np.where(np.isnan(a))])
print(a[np.where(a>12)])
a2    = np.array([1,2,3,4,5])
inds = np.where(a2>2)
print(a2[inds])
b = np.array([[[1, 2, 3],[4, 5, 6]],
            [[7, 8, 9],[10, 11, 12]]]) #3D matrix
print(b[1 ,0 ,2])
