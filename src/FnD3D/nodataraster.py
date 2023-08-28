# -*- coding: utf-8 -*-
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


Function to filter the no data raster
"""
import os
from osgeo import gdal

def nodataraster(prep_out):
    dss = gdal.Open(prep_out)
    bands = dss.GetRasterBand(1)
    statss = bands.GetStatistics(True,True)
    if statss[0]!=statss[1]:
        del(dss)
        del(bands)
        del(statss)
        del(prep_out)
    else:    
        del(dss)
        del(bands)
        del(statss)
        prep_del = prep_out
        del(prep_out)
        os.remove(str(prep_del))

