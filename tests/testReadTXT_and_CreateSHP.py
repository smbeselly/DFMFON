# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:49:26 2021

@author: sbe002

This script will serve the task 2.7

"""
import os
import pandas as pd
import glob

#TODO

# Read output file from each run folder

# compile the read output (txt) and concatenate

# Convert to txt back with header or as dataframe only


#source: https://stackoverflow.com/questions/20906474/import-multiple-csv-files-into-pandas-and-concatenate-into-one-dataframe
filepaths = [f for f in os.listdir(".") if f.endswith('.csv')]
df = pd.concat(map(pd.read_csv, filepaths))

# df = pd.DataFrame({'id': id_id,
#                        'speciesName' : speciesName}) add header