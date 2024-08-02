# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 23:23:00 2021

@author: sbe002

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

