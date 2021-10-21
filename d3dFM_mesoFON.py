# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Script to Coupling Delft3D and MesoFON

# Import the necessary packages

## Initiate the BMI Delft3D

## Read the D3D Domain and Map it to create the 'world' for MesoFON

## Get the limit of the approximated LLWL and build the outer line
# This is as the base for 'salinity approach' to kill the trees that are 
# located outside the allowed area

# =============================================================================
# # ~ try splitting the shp files and the raster
# # ~ create function based on this recipe "splitD3Dmap"
# Already done in testSplitD3D
# It creates a tiles of raster file with specific interval in meters or pixel
# =============================================================================


### Loop

    ## Read the trees location and create the polygon for Trachytope
    
    ## Do the first run of Delft3D
    
    # Pause the simulation to a certain delta t coupling (3 months)
    
    ## Script to process the output with BMI? or other default Python processing
    ## package for DFM
    
    ## Retrieve the D3DFM map based on WoO
    # run the "splitD3Dmap"
        # It will clip the map based on the master shapefile
        # Create raster based on the WoO as a Salinity, on the masterClipped 
        # However, it will only mark the unhabitable area as extreme salinity
        # ----
        # Split the WoO Map
        
    # the output will be the tiles of raster files with extreme value of salinity 
    # that limits the growth of the existing trees and kill the sapling directly
    
    ## Call the MesoFON and run with the newly created raster world
    # pause after 1 step (=3 months)
    # ~ MesoFON will pause and wait for the new raster files from Delft3D
    
    # receive the new trees location signal from MesoFON
    
    ## Rebuild the trees based on the tiles 

### End Loop

## Post Processing for the Delft3D Simulation
## Post Processing for the MesoFON Simulation  


