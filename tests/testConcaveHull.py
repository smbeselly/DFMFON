# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 15:49:09 2021

@author: sbe002
"""
# how to import the written function
# https://stackoverflow.com/questions/33409424/import-own-py-files-in-anaconda-spyder

import sys
print(sys.path)
sys.path.append('FND3D')
import numpy as np
import ConcaveHull as CH
############################################################
## tests
############################################################

nx = 5
ny = nx
tpoints = []
tpointsy = []
x = np.arange(nx)
for i in x:
        tpoints.append((i,0))
        tpointsy.append((0,i))
tpoints = np.asarray(tpoints)
tpointsy = np.asarray(tpointsy)

def test_GetNearestNeighbors_1d_x_00():
    neighbors = CH.GetNearestNeighbors(tpoints, (0,0),nx)
    assert np.array_equal(neighbors,tpoints)
def test_GetNearestNeighbors_1d_x_nx():
    neighbors = CH.GetNearestNeighbors(tpoints, (nx,nx),nx)
    assert np.array_equal(neighbors,tpoints[::-1])

def test_GetNearestNeighbors_1d_y_00():
    neighbors = CH.GetNearestNeighbors(tpointsy, (0,0),nx)
    assert np.array_equal(neighbors,tpointsy)
def test_GetNearestNeighbors_1d_y_nx():
    neighbors = CH.GetNearestNeighbors(tpointsy, (nx,nx),nx)
    assert np.array_equal(neighbors,tpointsy[::-1])


################################################
clock = np.array([[3,2], [4,2], [4,3], [4,4], [3,4], [2,4], [2,3], [2,2]])

def check_SortByAngle(points,currentPoint,prevPoint,i):
    sortedPoints = CH.SortByAngle(points,currentPoint,prevPoint)
    #for j in np.arange(len(clock)):
    #    print sortedPoints[j], clock[j], np.roll(clock,-i,axis=0)[j]
    assert np.array_equal(sortedPoints, np.roll(clock,-i,axis=0))

def test_SortByAngle_clock():
    i=0
    for point in clock:
        yield check_SortByAngle, clock, (3,3), point, i
        i=i+1

################################################
# simple test dataset
points = np.array([[10,  9], [ 9, 18], [16, 13], [11, 15], [12, 14], [18, 12],
                   [ 2, 14], [ 6, 18], [ 9,  9], [10,  8], [ 6, 17], [ 5,  3],
                   [13, 19], [ 3, 18], [ 8, 17], [ 9,  7], [ 3,  0], [13, 18],
                   [15,  4], [13, 16]])
points_solution_k_5 = np.array([[3, 0],[10,  8],[15,  4],[18, 12],[13, 18],[13, 19],
                               [ 9, 18],[6, 18],[3, 18],[2, 14],[9, 9],[5, 3],[3, 0]
                               ])
hull = CH.concaveHull(points, 5)
CH.plotPath(points, hull)

# points to test what happens if all points intersect
points_intersect = np.array([[1,1],[10,3],[11,8],[9,14],[15,21],[-5,15],[-3,10],
                            [2,5],    # from here the distracting points
                            [9,10],[8,9],[8,11],[8,12],[9,11],[9,12]
                            ])
points_intersect_solution = np.array([[1, 1],[10,  3],[11,  8],[9, 14],[15, 21],
                                     [-5, 15],[-3, 10],[1, 1]
                                     ])

hull2 = CH.concaveHull(points_intersect, 5)
CH.plotPath(points_intersect, hull2)


points_E = np.array([[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[6,2],[6,3],[5,3],[4,3],
                    [3,3],[3,4],[3,5],[4,5],[5,5],[5,6],[5,7],[4,7],[3,7],[3,8],
                    [3,9],[4,9],[5,9],[6,9],[6,10],[6,11],[5,11],[4,11],[3,11],[2,11],
                    [1,11],[1,10],[1,9],[1,8],[1,7],[1,6],[1,5],[1,4],[1,3],[1,2],
                    [5,2],[4,2],[3,2],[2,2],[2,3],[2,4],[2,5],[2,6],[2,7],[2,8],
                    [2,9],[2,10],[3,10],[4,10],[5,10],[3,6],[4,6],[5,6],[4.5,7],[3,8.5],
                    ])

hull3 = CH.concaveHull(points_E, 3)
CH.plotPath(points_E, hull3)