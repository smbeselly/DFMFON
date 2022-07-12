## Test Vegetation with vege_corner.pli

import numpy as np
import bmi.wrapper
import os
from dfm_tools.io.polygon import Polygon
import matplotlib.path as mpltPath

# define modeling environment
D3D_HOME = os.path.join(r'D:\Git\d3d_meso\Model-Execute\D3DFM\2sebrian_20220518')
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
dflowfm_path = os.path.join(D3D_HOME, 'dflowfm','bin')
dflowfm_engine = os.path.join(dflowfm_path, 'dflowfm.dll')

D3D_workdir = os.path.join(r'D:\Git\d3d_meso\Model-Execute\D3DFM\model_run_6_flat')
config_file = os.path.join(D3D_workdir, 'model_run_6_flat.xml') # funnel morpho
mdu_file = os.path.join(D3D_workdir,'dflowfm', 'FlowFM.mdu')

# D3D_workdir = os.path.join(r'D:\Git\d3d_meso\Model-Execute\D3DFM\model_run_6_flat')
# config_file = os.path.join(D3D_workdir, 'model_run_6_flat.xml') # funnel morpho
# mdu_file = os.path.join(D3D_workdir,'dflowfm', 'FlowFM.mdu')

### Initialize BMI Model
## Add corrects locations to environment variable PATH for DFM
os.environ['PATH'] = os.path.join(D3D_HOME, 'share', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dflowfm', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dimr', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dwaves', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'esmf', 'scripts') \
+ ";" + os.path.join(D3D_HOME, 'swan', 'scripts')

## Define DFM wrapper
# os.chdir(dflowfm_path)
model_dfm = bmi.wrapper.BMIWrapper(engine=dflowfm_engine, configfile=mdu_file)

## Define and initialise DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()

#%% test
## initialisation of vegetation variables
rnveg=model_dfm.get_var('rnveg')  # dichtheid
diaveg=model_dfm.get_var('diaveg') # diameter
stemheight=model_dfm.get_var('stemheight') # hoogte
ndx = model_dfm.get_var_shape('s1')

Cdbase=5
drag_coeff = np.random.random()*Cdbase*np.ones(ndx)
coupling_period_model = 43200

# Read shear polygon, and get indices
xz=model_dfm.get_var('xz')
yz=model_dfm.get_var('yz')
# pli = Polygon.fromfile(os.path.join(D3D_workdir,'dflowfm', 'vege.pli'))
# pli = Polygon.fromfile(os.path.join(D3D_workdir,'dflowfm', 'vege_corner.pli'))
pli = Polygon.fromfile(os.path.join(D3D_workdir,'dflowfm', 'vege_tim_lauttt.pli'))

path = mpltPath.Path(pli[0][0])
ind = path.contains_points(np.transpose((xz,yz))).astype(int)

rnveg=np.full(ndx,1.0)
diaveg=np.full(ndx,0.2)
stemheight=np.full(ndx,7.0)

model_dfm.set_var('rnveg',ind*rnveg)
model_dfm.set_var('diaveg',ind*diaveg)
model_dfm.set_var('stemheight',ind*stemheight)

#%% run the model
t=0

    # update the variable with new value
cdf=np.random.random()*Cdbase
drag_coeff = np.ones(ndx)*cdf*ind 
model_dfm.set_var('Cdvegsp',drag_coeff)
print('drag_coeff set to ',cdf)
while t<coupling_period_model:
    # model_dfm.set_var('s1',water_level)
    # cdf=np.random.random()*Cdbase
    # drag_coeff = np.ones(ndx)*cdf*ind   
    # model_dfm.set_var('Cdvegsp',drag_coeff)
    # print('drag_coeff set to ',cdf)
    model_dfm.update()
    dts = model_dfm.get_time_step()
    t=t+dts
    
#Finalize the running
model_dfm.finalize() 
