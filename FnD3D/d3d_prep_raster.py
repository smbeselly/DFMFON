# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:12:08 2021

@author: sbe002

This is the accompaying function to process the preparation for D3DFM-MesoFON
"""

import os
import numpy as np
import numpy.ma as ma
from dfm_tools.get_nc import get_netdata
from shapely.geometry import mapping, Polygon
from osgeo import gdal
import ConcaveHull as CH

#%% Create ConcaveHull for nc_input

def d3dConcaveHull(nc_in,k_hull):
    """Read the UGrid = mesh2d_node_x and mesh2d+node_y
    and create ConcaveHull Polygon
    
    Parameters
    ---------
    nc_in = netCDF UGrid file (both old and new UGrid)
    k_hull = k value to build Concave Hull
    
    Returns
    ---------
    Concave Hull Polygon Object
    """

    file_nc = nc_in
    
    ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
    matrix = np.block([[ma.compressed(ugrid_all.mesh2d_node_x)],[ma.compressed(ugrid_all.mesh2d_node_y)],[ma.compressed(ugrid_all.mesh2d_node_z)]]).T
    
    matrix= (np.delete(matrix, np.where(matrix == -999)[0], axis=0))
    matrix = (matrix[~np.isnan(matrix).any(axis=1)])
    
    mat_hull = np.block([[matrix[:,0]],[matrix[:,1]]]).T # create matrix x,y for hull
    
    #% creating concave hull
    hull = CH.concaveHull(mat_hull, k_hull) # utk case ini k=9
    # CH.plotPath(matrix, hull)
    
    # create polygon
    poly = Polygon(hull)
    # plt.plot(*poly.exterior.xy)
    
    return poly

#%% Create SHP from Polygon
def d3dPolySHP(poly_in, dir_out, out_poly_name, projection):
    """Create Shapefile from shapely polygon feature
    
    Parameters
    ---------
    poly_in = Polygon file
    dir_out = Output directory
    out_poly_name = name of the shapefile as string
    projection = EPSG Code as string
    
    Returns
    ---------
    Shapefile with 'out_poly_name.shp'
    """
    import fiona
    from fiona.crs import from_epsg
    # output location
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    schema = {
        'geometry': 'Polygon',
        'properties': {'id': 'int', 'area':'float',
                       'length':'float'},
    }
    # Write a new Shapefile
    # source: https://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
    with fiona.open((str(dir_out)+str('\\')+str(out_poly_name)+'.shp'), 'w', \
                    'ESRI Shapefile', crs=from_epsg(projection), schema=schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        c.write({
            'geometry': mapping(poly_in),
            'properties': {'id': 123, 'area':poly_in.area,
                           'length':poly_in.length},
        })

#%% Create CSV build raster and clip
def d3dCSV2ClippedRaster(concave_path, concave_name, EPSG_coord, matrix, x_res, y_res, no_data_val, shp_clip, affix):
    """Create raster from the np.array/ matrix and clip the raster based on shapefile
        
    Parameters
    ---------
    concave_path = path to the output folder
    concave_name = name of the created raster in string
    EPSG_coord = EPSG Code as string
    matrix = data matrix as np.array in x,y,z
    x_res = x axis resolution in meter
    y_res = y axis resolution in meter
    no_data_val = define no data value
    shp_clip = refer to the shapefile source as mask
    affix = define the affix of the clipped raster in string
    
    Returns
    ---------
    Raster from Matrix and
    Clipped Raster
    """

    np.savetxt(str(concave_path)+str('\\')+str(concave_name)+'.csv', matrix, fmt="%f", delimiter=",", header='Lon,Lat,Alt',comments='')
    
    # Create .vrt file based on the created .csv
    lines = ['<OGRVRTDataSource>', '<OGRVRTLayer name='+'"'+str(concave_name)+'"'+'>',
             '<SrcDataSource>'+str(concave_name)+".csv"+'</SrcDataSource>',
             '<GeometryType>wkbPoint</GeometryType>',
             '<LayerSRS>EPSG:'+str(EPSG_coord)+'</LayerSRS>',
             '<GeometryField separator=" " encoding="PointFromColumns" x="Lon" y="Lat" z="Alt"/>',
             '</OGRVRTLayer>',
             '</OGRVRTDataSource>']
    with open(str(concave_path)+str('\\')+str(concave_name)+'.vrt', 'w') as f:
        f.write('\n'.join(lines))
        
    # Call and run gdal_grid to create raster file from bathymetry
    # abs_path = os.path.join('D:/Git/d3d_meso/')
    os.chdir(concave_path)
    vrt_in = str(concave_name)+'.vrt'
    csv_in = str(concave_name)
    ras_out = str(concave_name)+'.tif'
    x_min = np.min(matrix[:,0])
    x_max = np.max(matrix[:,0])
    y_min = np.min(matrix[:,1])
    y_max = np.max(matrix[:,1])
    
    # x_res = 10
    # y_res = 10
    
    # no_data_val = -999
    
    # create a raster with linear interpolation     
    command = 'gdal_grid -zfield "Alt" -a linear:radius=0.0:nodata={no_data_val} -ot Float32 \
                -txe {x_min} {x_max} -tye {y_min} {y_max} \
                -tr {x_res} {y_res} -l {csv_in} {vrt_in} {ras_out} \
                --config GDAL_NUM_THREADS ALL_CPUS'
    
    os.system(command.format(no_data_val=no_data_val, x_min=x_min, x_max=x_max, 
                             y_min=y_min, y_max=y_max, x_res=x_res, y_res=y_res, 
                             csv_in=csv_in, vrt_in=vrt_in, ras_out=ras_out))
    
    # call and run gdalwarp to clip the raster and add no data value
    cut_cl = str(shp_clip)
    cut_call = str(shp_clip)+'.shp'
    ras_clip = str(concave_name)+affix+'.tif'
    
    command_warp = 'gdalwarp -overwrite -of GTiff -cutline {cut_call} -cl {cut_cl} \
                    -crop_to_cutline -dstnodata {no_data_val} {ras_out} {ras_clip}'
    os.system(command_warp.format(cut_call=cut_call, cut_cl=cut_cl, no_data_val=no_data_val, 
                                  ras_out=ras_out, ras_clip=ras_clip))

#%% Tile the raster and filter + delete nodata tiled raster
 
def d3dRaster2Tiles(out_path, output_filename, ras_clip, tile_size_x, tile_size_y, CreateSHP=True):
    """Create tiles of raster tiles of shapefile
        
    Parameters
    ---------
    out_path = path to save the file
    output_filename = name of the tiled rasters in string
    ras_clip = original clipped raster to be tiled
    tile_size_x = resolution of the tile in x direction in meter
    tile_size_y = resolution of the tile in y direction in meter
    CreateSHP = boolean, set True if you want to create SHP, set False for no SHP
    
    
    Returns
    ---------
    Tiles of the clipped raster
    and tiles of shapefile of the clipped raster if CreateSHP=True
    """
    ds = gdal.Open(ras_clip)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    # stats = band.GetStatistics(True,True) # results = min, max, mean, StdDev
    # stats[1]-stats[0] = max-min
    xsize = band.XSize
    ysize = band.YSize
    # get coordinates of upper left corner
    # xmin = gt[0]
    # ymax = gt[3]
    res = gt[1]
    
    # determine total length of raster
    # xlen = res * ds.RasterXSize # convert the size in meter to number of pixels
    # ylen = res * ds.RasterYSize
    
    # size of a single tile
    xsize_tile = round(tile_size_x/res) #num of pixels in tile_size_x
    ysize_tile = round(tile_size_y/res) ##num of pixels in tile_size_y
    
    # Tile the raster domain as the prefered tile size in meter
    for i in range(0, xsize, xsize_tile):
        for j in range(0, ysize, ysize_tile):
            prep_out = str(out_path) +str('\\')+ str(output_filename) + str(i) + "_" + str(j) + ".tif"        
            command = 'gdal_translate -of GTIFF -srcwin {i}, {j}, {xsize_tile}, \
                {ysize_tile} {prep_conv_out} {prep_out}'
            os.system(command.format(i=i, j=j, xsize_tile=xsize_tile, ysize_tile=ysize_tile, 
                                     prep_conv_out=ras_clip, prep_out=prep_out))
            # below is to filter and delete the grid with all nodata value
            dss = gdal.Open(prep_out)
            bands = dss.GetRasterBand(1)
            statss = bands.GetStatistics(True,True)
            
            if statss[0] == 0 and statss[1] == 0:
                del(dss,bands,statss)
                prep_del = prep_out
                del(prep_out)
                os.remove(str(prep_del))
            else:
                #build shapefile
                if CreateSHP == True:
                    import rasterio as rio
                    dt = rio.open(prep_out)
                    crs_raster= dt.crs
                    left,bottom, right, top = dt.bounds #[left, bottom, right, top]
                    # Create Shapefile  
                    from shapely.geometry import box
                    poly_ras= box(left, bottom, right, top) #create poly box
                    # plt.plot(*poly_ras.exterior.xy)
                    out_poly_namet = str(output_filename) + str(i) + "_" + str(j)
                    d3dPolySHP(poly_ras, out_path, out_poly_namet, str(crs_raster)[5:])
                    del(dt, crs_raster, left, bottom, right, top, poly_ras)
                    # import gc 
                    # gc.collect()
                # deleting variables
                del(dss,bands,statss,prep_out)
            
            # if statss[0]!=statss[1]:
            #     #build shapefile
            #     if CreateSHP == True:
            #         import rasterio as rio
            #         dt = rio.open(prep_out)
            #         crs_raster= dt.crs
            #         left,bottom, right, top = dt.bounds #[left, bottom, right, top]
            #         # Create Shapefile  
            #         from shapely.geometry import box
            #         poly_ras= box(left, bottom, right, top) #create poly box
            #         # plt.plot(*poly_ras.exterior.xy)
            #         out_poly_namet = str(output_filename) + str(i) + "_" + str(j)
            #         d3dPolySHP(poly_ras, out_path, out_poly_namet, str(crs_raster)[5:])
            #         del(dt, crs_raster, left, bottom, right, top, poly_ras)
            #         # import gc 
            #         # gc.collect()
            #     # deleting variables
            #     del(dss,bands,statss,prep_out)
			
			# else:    
            #     del(dss,bands,statss)
            #     prep_del = prep_out
            #     del(prep_out)
            #     os.remove(str(prep_del))

#%% Tile the raster and filter + retain nodata tiled raster               
def d3dRaster2TilesRetain(out_path, output_filename, ras_clip, tile_size_x, tile_size_y, CreateSHP=True):
    """Create tiles of raster tiles of shapefile
        
    Parameters
    ---------
    out_path = path to save the file
    output_filename = name of the tiled rasters in string
    ras_clip = original clipped raster to be tiled
    tile_size_x = resolution of the tile in x direction in meter
    tile_size_y = resolution of the tile in y direction in meter
    CreateSHP = boolean, set True if you want to create SHP, set False for no SHP
    
    
    Returns
    ---------
    Tiles of the clipped raster
    and tiles of shapefile of the clipped raster if CreateSHP=True
    """
    ds = gdal.Open(ras_clip)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    # stats = band.GetStatistics(True,True) # results = min, max, mean, StdDev
    # stats[1]-stats[0] = max-min
    xsize = band.XSize
    ysize = band.YSize
    # get coordinates of upper left corner
    # xmin = gt[0]
    # ymax = gt[3]
    res = gt[1]
    
    # determine total length of raster
    # xlen = res * ds.RasterXSize # convert the size in meter to number of pixels
    # ylen = res * ds.RasterYSize
    
    # size of a single tile
    xsize_tile = round(tile_size_x/res) #num of pixels in tile_size_x
    ysize_tile = round(tile_size_y/res) ##num of pixels in tile_size_y
    
    # Tile the raster domain as the prefered tile size in meter
    for i in range(0, xsize, xsize_tile):
        for j in range(0, ysize, ysize_tile):
            prep_out = str(out_path) +str('\\')+ str(output_filename) + str(i) + "_" + str(j) + ".tif"        
            command = 'gdal_translate -of GTIFF -srcwin {i}, {j}, {xsize_tile}, \
                {ysize_tile} {prep_conv_out} {prep_out}'
            os.system(command.format(i=i, j=j, xsize_tile=xsize_tile, ysize_tile=ysize_tile, 
                                     prep_conv_out=ras_clip, prep_out=prep_out))
            # below is to filter and delete the grid with all nodata value
            dss = gdal.Open(prep_out)
            bands = dss.GetRasterBand(1)
            statss = bands.GetStatistics(True,True)
            # if statss[0]!=statss[1]:
            #build shapefile
            if CreateSHP == True:
                import rasterio as rio
                dt = rio.open(prep_out)
                crs_raster= dt.crs
                left,bottom, right, top = dt.bounds #[left, bottom, right, top]
                # Create Shapefile  
                from shapely.geometry import box
                poly_ras= box(left, bottom, right, top) #create poly box
                # plt.plot(*poly_ras.exterior.xy)
                out_poly_namet = str(output_filename) + str(i) + "_" + str(j)
                d3dPolySHP(poly_ras, out_path, out_poly_namet, str(crs_raster)[5:])
                del(dt, crs_raster, left, bottom, right, top, poly_ras)
                # import gc 
                # gc.collect()
            # deleting variables
            del(dss,bands,statss,prep_out)					   																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																					
            
