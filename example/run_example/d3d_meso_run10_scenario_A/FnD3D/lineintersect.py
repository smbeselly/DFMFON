#!/bin/env python
"""
Copyright notice
-------------------------
This library is developed as part of the PhD research of
Sebrian Mirdeklis Beselly Putra conducted at IHE Delft Institute 
for Water Education and Delft University of Technology

The author  of this library is:
    Sebrian Beselly
    s.besellyputra@un-ihe.org
    s.m.beselly@tudelft.nl
    sebrian@ub.ac.id
    
    IHE Delft Institute for Water Education,
    PO Box 3015, 2601DA Delft
    the Netherlands
    
This library is free software: you can redistribute it and/or modify 
it under the terms of the GPL-3.0 license
    
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GPL-3.0 license for more details.

Publication related to this library
Beselly, S.M., U. Grueters, M. van Der Wegen, J. Reyns, J. Dijkstra, and D. Roelvink. “Modelling Mangrove-Mudflat Dynamics with a Coupled Individual-Based-Hydro-Morphodynamic Model.” Environmental Modelling & Software, August 28, 2023, 105814. https://doi.org/10.1016/j.envsoft.2023.105814.


"""
##
## determine if two line segments intersect
## see: martin-thoma.com/how-to-check-if-two-line-segments-intersect/
##

import numpy as np

def doBoundingBoxesIntersect(a, b, c, d):
    '''
    Check if bounding boxes do intersect. If one bounding box touches
    the other, they do intersect.
    First segment is of points a and b, second of c and d.
    '''
    ll1_x = min(a[0],b[0]); ll2_x = min(c[0],d[0])
    ll1_y = min(a[1],b[1]); ll2_y = min(c[1],d[1])
    ur1_x = max(a[0],b[0]); ur2_x = max(c[0],d[0])
    ur1_y = max(a[1],b[1]); ur2_y = max(c[1],d[1])

    return ll1_x <= ur2_x and \
           ur1_x >= ll2_x and \
           ll1_y <= ur2_y and \
           ur1_y >= ll2_y

def isPointOnLine(a,b,c):
    '''
    Check if a point is on a line.
    '''
    # move to origin
    aTmp = (0,0)
    bTmp = (b[0] - a[0], b[1] - a[1])
    cTmp = (c[0] - a[0], c[1] - a[1])
    r = np.cross(bTmp, cTmp)
    return np.abs(r) < 0.0000000001

def isPointRightOfLine(a,b,c):
    '''
    Check if a point (c) is right of a line (a-b).
    If (c) is on the line, it is not right it.
    '''
    # move to origin
    aTmp = (0,0)
    bTmp = (b[0] - a[0], b[1] - a[1])
    cTmp = (c[0] - a[0], c[1] - a[1])
    return np.cross(bTmp, cTmp) < 0

def lineSegmentTouchesOrCrossesLine(a,b,c,d):
    '''
    Check if line segment (a-b) touches or crosses
    line segment (c-d).
    '''
    return isPointOnLine(a,b,c) or \
           isPointOnLine(a,b,d) or \
          (isPointRightOfLine(a,b,c) ^
           isPointRightOfLine(a,b,d))

def doLinesIntersect(a,b,c,d):
    '''
    Check if line segments (a-b) and (c-d) intersect.
    '''
    return doBoundingBoxesIntersect(a,b,c,d) and \
           lineSegmentTouchesOrCrossesLine(a,b,c,d) and \
           lineSegmentTouchesOrCrossesLine(c,d,a,b)


##############################
## Tests
##############################

def test_doBoundingBoxesIntersect():
    A=(1,1); B=(2,2); C=(3,1); D=(4,2)
    assert doBoundingBoxesIntersect(A,B,C,D) == False
    A=(1,2); B=(3,4); C=(2,1); D=(4,3)
    assert doBoundingBoxesIntersect(A,B,C,D) == True

def test_isPointOnLine():
    A=(1,1); B=(3,3); C=(2,2)
    assert isPointOnLine(A,B,C) == True
    A=(1,1); B=(3,3); C=(3,2)
    assert isPointOnLine(A,B,C) == False


def test_isPointRightOfLine():
    A=(1,1); B=(3,3); C=(2,2)
    assert isPointRightOfLine(A,B,C) == False
    A=(1,1); B=(3,3); C=(3,2)
    assert isPointRightOfLine(A,B,C) == True
    A=(1,1); B=(3,3); C=(1,2)
    assert isPointRightOfLine(A,B,C) == False

# a lot of testcases to be tested with the final function

def tcase(name):
    if name == 'F1':
        return (0,0), (7,7), (3,4), (4,5), False
    elif name == 'F2':
        return (-4,4), (-2,1), (-2,3), (0,0), False
    elif name == 'F3':
        return (0,0), (0,1), (2,2), (2,3), False
    elif name == 'F4':
        return (0,0), (0,1), (2,2), (3,2), False
    elif name == 'F5':
        return (-1,-1), (2,2), (3,3), (5,5), False
    elif name == 'F6':
        return (0,0), (1,1), (2,0), (0.5,2), False
    elif name == 'F7':
        return (1,1), (4,1), (2,2), (3,2), False
    elif name == 'F8':
        return (0,5), (6,0), (2,1), (2,2), False
    elif name == 'T1':
        return (0,-2), (0,2), (-2,0), (2,0), True
    elif name == 'T2':
        return (5,5), (0,0), (1,1), (8,2), True
    elif name == 'T3':
        return (-1,0), (0,0), (-1,-1), (-1,1), True
    elif name == 'T4':
        return (0,2), (2,2), (2,0), (2,4), True
    elif name == 'T5':
        return (0,0), (5,5), (1,1), (3,3), True
    elif name == 'T6':
        return (0,0), (3,3), (0,0), (3,3), True

cases = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8',
         'T1', 'T2', 'T3', 'T4', 'T5', 'T6']

def check_intersection(name):
    A,B,C,D, result = tcase(name)
    assert doLinesIntersect(A,B,C,D) == result

def test_doLinesIntersect():
    for case in cases:
        yield check_intersection, case

