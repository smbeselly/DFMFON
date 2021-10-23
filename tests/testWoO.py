# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 16:26:23 2021

@author: sbe002

still remain the same as the old script
"""

## This script is based on the Broek_script.py
# Location: D:\BOOKS\IHE_Delft\MSc Thesis-Other\Broek_script.py
import numpy as np
import bmi.wrapper
import os
import datetime
from scipy import integrate
import faulthandler
faulthandler.enable()

## Specify paths of dll-files and input-files
D3D_HOME = r'E:\Data Delft3Dmodel\Deltversion\Marconiv2\MarconiModel\MarconiModel\Marconi_v02\Code1709\windows\oss_artifacts_x64_63721\x64'
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')

dflowfm_path = os.path.join(D3D_HOME, 'dflowfm', 'bin', 'dflowfm.dll')
workdir = r"E:\Data Delft3Dmodel\Deltversion\Marconiv2\MarconiModel\MarconiModel\Marconi_v02"
config_file = os.path.join(workdir, 'dimr_config.xml')
mdu_file = os.path.join(workdir, 'fm', 'Marconi_v02.mdu')
grid_file = os.path.join(workdir, 'fm', 'grid_marconi_v07_net.nc')
figsavefolder= r"E:\Data Delft3Dmodel\Deltversion\Marconiv2\MarconiModel\MarconiModel\Marconi_v02\figs"

print("D3D_HOME :" + D3D_HOME)
print("dimr_path :" + dimr_path)
print("config_file:" + config_file)

## Add corrects locations to environment variable PATH
os.environ['PATH'] = os.path.join(D3D_HOME, 'share', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dflowfm', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dimr', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'dwaves', 'bin') \
+ ";" + os.path.join(D3D_HOME, 'esmf', 'scripts') \
+ ";" + os.path.join(D3D_HOME, 'swan', 'scripts')

## Define DFM wrapper
model_dfm = bmi.wrapper.BMIWrapper(engine=dflowfm_path, configfile=mdu_file)

## Define and initialise DIMR wrapper
model_dimr = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
model_dimr.initialize()
print ('model initialized')

## Get the pointers to important model variables of FlowFM, list of accessible parameters and meaning can be found at:https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/dflowfm/packages/dflowfm_lib/include/bmi_get_var.inc

ndx=model_dfm.get_var('ndx')
ndxi=model_dfm.get_var('ndxi')
xzw=model_dfm.get_var('xzw')
yzw=model_dfm.get_var('yzw')
lnx=model_dfm.get_var('lnx')
lnxi=model_dfm.get_var('lnxi')
ln=model_dfm.get_var('ln')
dx=model_dfm.get_var('dx')
wu=model_dfm.get_var('wu')
ba=model_dfm.get_var('ba')
bedlevel=model_dfm.get_var('bl')

## initialisation of vegetation variables
rnveg=model_dfm.get_var('rnveg')
stemdia=model_dfm.get_var('diaveg')
stemheight=model_dfm.get_var('stemheight')

## settings
# define number of timesteps (nt) and duration of time steps (tstep) in [s]
mtpervt=12*3600 #model run time (s) between vegetation updates/checks
daysecs=3600*24
MF=4.0
tstep_season = int(180/(mtpervt*MF/(daysecs))) #43 # length of a season is 180 days / 4.16666 days = 43 days
nt=tstep_season #or 1 for a test
model_days = nt*mtpervt / (3600*24)
model_months = nt*mtpervt / (86400*30)
ncells = len(rnveg)
print ('models spans ', model_days, ' days, or ', model_months, ' months')
print ('# of vegetation timesteps ', nt)
print ('timestep length ', mtpervt ,' seconds')

## WoO parameters
WoO_T = np.zeros(shape=(len(bedlevel),)) # WoO progression tracking array: how far is each grid cells into the WoO
W2steps = np.zeros(shape=(len(bedlevel),)) # How many timesteps in second window completed
hcrit = 0.01 # critical inundation depth [m]
CED1 = 8e-03 # initial critical erosion depth single timestep [m]
CED2 = 27e-03 # mature critical erosion depth single timestep [m]
Savg = 15*((1.0)/(MF*7*24))*1.0e-03 #originally /by morfac!# critical average sedimentation over plant life in [m/s] (15 mm / week)
Smin = 1*60*60 # minimum age for average sedimentation chekc in [s]
Eavg = 5*((1.0)/(MF*7*24))*1.0e-03 #originally /by morfac!# critical average erosion over plant life [m/s] (5 mm / week)
Emin = 1*60*60 # minimum age for average erosion check in [s]
alpha = 1.50
height_W2=0.02
PW2=20 #Stem density at the start of W2

import array as arr
EP = arr.array('d', [0.0075, 0.0075, 0.0075, 0.0075, 0.0075, 0.0075]) #Schwarz: 0.92, 0.92, 0.84, 0.76, 0.76, 0.76
W1 = (2.5*daysecs)/MF #Divided by morfac length of Window 1 in [s] (2.5 days default)
wdcrit = 0.01 # critical inundation depth [m]
W2 = (90*daysecs)/MF #Divided by morfac length of Window 2 in [s]; days*daylength (90 days default)

## preallocating the matrix used for saving historic bed level values and vegetation density values
blh = np.zeros(shape=(len(bedlevel),3))
rnveg_historic = np.zeros(shape=(10553,((nt))))
bl_historic = np.zeros(shape=(len(bedlevel),((nt))))
stemheight_historic=np.zeros(shape=(len(bedlevel),((nt))))
EPSuc=np.zeros(shape=(len(bedlevel)))
Veghist = np.zeros(shape=(len(xzw),1))
Plantonce=np.zeros(shape=(len(bedlevel)))
c=np.zeros(shape=(len(bedlevel)))

## population dynamics parameters
yearsec=24*60*60*365
P0=60.0 # Plant density of seedlings after W2 (stems/m2)
r=20.0/(365*2) # intrinsic growth rate of plant density (per year)
K=600.0 # carrying capacity of plant density (stems/m2)
Ctau=(30.0/(365*24))*yearsec # plant erosion coefficient due to bed shear stress (m^-2 s^1)
tau_crp=0.25 # critical bed shear stress for plant erosion (N/m2)
Cinund=3000/(365*24) # plant erosion coefficient due to inundation stress (m^-3 y^-1)
inund_crp=1.0 # critical inundation height at high tide for plant erosion (m)
diam_def=0.005 # default diameter stems
height_def=0.28 # default height stems
D=0 # diffusion coefficient [m^2 y^-1]

stemdia[:]=diam_def
model_dfm.set_var('diaveg',stemdia)
print ('stemdia set: ',stemdia)
stemheight[:]=height_def
model_dfm.set_var('stemheight',stemheight)
print ('stemheight set: ',stemheight)
Vegsectionsnodes=np.loadtxt("E:\Data Delft3Dmodel\Deltversion\Marconiv2\MarconiModel\MarconiModel\Marconi_v02\CellsSectionsv2.txt",dtype=int) #The sown seeds at Sections E, F and G
rnveg[:]=0
model_dfm.set_var('rnveg',rnveg)
print ('rnveg set: ' + str(rnveg))

## locating indices of possible vegetated gridcells used for the WOO loop (it is unnecessary to check for windows in grid cell located below -1m, no vegetation can grow there)
count=0
veg_area = np.zeros(shape=(10553,))
for j in range(0,10553):
	if bedlevel[j]>=-1 and bedlevel[j]<3: #NB not exact numbers, just an approximation
		veg_area[count]=int(j)
		count=count+1
veg_area=veg_area[:count]

def makeinirandfield (P):
	P[range(ndxi)]=0
	return(P)
	
def erosveg(P,tau_max):
	dptau=-P*Ctau*(tau_max-tau_crp)
	dptau[dptau>0.0]=0.0
	return(dptau)
	
def inundveg(P,wd_mn):
	dpinund=-P*Cinund*(wd_mn-inund_crp)
	dpinund[dpinund>0.0]=0.0
	return(dpinund)
	
def vegRoC (P,t0):
	dPdt=np.zeros(ndxi)
	# logistic growth
	dPdt=dPdt+r*P*(1-P/K)
	# Diffusion
	dPdt=dPdt+0
	# erostau
	dPdt=dPdt+erosveg(P,tau_mn)
	# inundation stress
	dPdt=(dPdt+inundveg(P,wd_mn))*MF
	# return rate of change
	return(dPdt)
	# population dynamics equation in winter season
def vegWinter (P,t0):
	dPdt=np.zeros(ndxi)
	# erostau
	dPdt=dPdt+erosveg(P,tau_mn)
	# inundation stress
	dPdt=dPdt+inundveg(P,wd_mn)
	# return rate of change
	dPdt=dPdt
	return(dPdt)
	
##set initial vegetation field
rnveg=makeinirandfield(rnveg[range(ndxi)]) # random seeding of plants
## model update and WoO loop
realtime=0
tdata=1
tt = 0
summersteps=0
wintersteps=0
year=0
t=np.array([0.0,1.0])

for i in range(1,nt+1):
	model_dimr.update(mtpervt) #update the model
	is_dtint=model_dfm.get_var('is_dtint') #timespan over which statistics are kept, needs to be reset after each vt (otherwise accumulative)
	print(is_dtint)
	realtime=realtime+is_dtint
	print('real time is ' + str(realtime) + ' s, or ' + str(realtime/daysecs) + ' days')
	print('step ' + str(i) + ' of ' + str(nt))
	bl=model_dfm.get_var('bl')
	hs=model_dfm.get_var('hs') #load in the water level
	blh[:,2]=bl #store the bed level in a matrix so it can be checked against previous bed levels in the WoO model
	
## summer time step calculations
	if summersteps < tstep_season:
		# population dynamic (growth,decay tau, decay wd) model every timestep
		is_maxvalsnd=model_dfm.get_var('is_maxvalsnd')
		is_sumvalsnd=model_dfm.get_var('is_sumvalsnd')
		bedlvl=model_dfm.get_var('bl')
		tau_mn=is_sumvalsnd[range(ndxi),0]/is_dtint
		vel_mn=is_sumvalsnd[range(ndxi),1]/is_dtint
		wd_mn=is_sumvalsnd[range(ndxi),2]/is_dtint
		tau_max=is_maxvalsnd[range(ndxi),0]
		vel_max=is_maxvalsnd[range(ndxi),1]
		wd_max=is_maxvalsnd[range(ndxi),2]
		
		# run plant growth model
		rslt=integrate.odeint(vegRoC, rnveg[range(ndxi)], t,h0=.05,hmax=.1)
		rslt=integrate.odeint(vegRoC, rnveg[range(ndxi)], t)
		for veg in range(ndxi):
			if WoO_T[veg]==-1:
				if rnveg[veg]>0.1: #decay due to 		bed shear stress
					rnveg[veg]=rslt[1,[veg]]
				if rnveg[veg]<=0.1 and WoO_T[veg]==-1 and rnveg[veg]!=0: #if established vegetation has vanished, WOO opens and is checked again
					rnveg[veg]=0.0
					WoO_T[veg]=0
				elif rnveg[veg]<=0.1 and rnveg[veg]!=0: #decay of diffused vegetation
					rnveg[veg]=0.0
		if year==0:
			for j in veg_area:
				k=int(j)
				c[k]=EP[year]
				a=np.random.rand(1)
				if Vegsectionsnodes[k]==1 and Plantonce[k]==0:
					Plantonce[k]=1
					c[k]=0.55;
				if a >=(1-c[k]) or EPSuc[k]==1:
					EPSuc[k]=1
					if WoO_T[k] == -1:
						continue
					elif WoO_T[k] <= W1:
						if hs[k] < hcrit:
							WoO_T[k]=WoO_T[k]+mtpervt 
						else:
							WoO_T[k]=0 
							EPSuc[k]=0
					elif WoO_T[k]> W1 and WoO_T[k] <= W1+W2:
						if W2steps[k]==0:
							blh[k,0]=blh[k,2]
							blh[k,1]=blh[k,2] 
							rnveg[k]=PW2 
							stemheight[k]=height_W2 
							W2steps[k]=1
						if (-blh[k,1] + blh[k,2]) > -(CED1 + ((WoO_T[k]-W1)/W2) * (CED2 - CED1) + (alpha * (-blh[k,0]+blh[k,2]))) and (-blh[k,0] + blh[k,2])/max((WoO_T[k]-W1),Emin) > -Eavg and (-blh[k,0] + blh[k,2])/max((WoO_T[k]-W1),Smin) < Savg:
							WoO_T[k]=WoO_T[k]+mtpervt 
							blh[k,1]=blh[k,2]
							stemheight[k]=height_W2 + (WoO_T[k]-W1) * (height_def-height_W2)/W2 
						else:
							WoO_T[k]=0 
							W2steps[k]=0 
							stemheight[k]=0 
							rnveg[k]=0 
							EPSuc[k]=0
					elif WoO_T[k] > W1+W2: 
						rnveg[k]=P0  
						WoO_T[k] = -1
						W2steps[k]=0
		elif year!=0:
			for j in veg_area:
				k=int(j) 
				a=np.random.rand(1)
				if a >=(1-(EP[year])) or EPSuc[k]==1:
					EPSuc[k]=1
					if WoO_T[k] == -1:
						continue
					elif WoO_T[k] <= W1:
						if hs[k] < hcrit:
							WoO_T[k]=WoO_T[k]+mtpervt 
						else:
							WoO_T[k]=0 
							EPSuc[k]=0
					elif WoO_T[k]> W1 and WoO_T[k] <= W1+W2: 
						if W2steps[k]==0:
							blh[k,0]=blh[k,2]
							blh[k,1]=blh[k,2]
							rnveg[k]=PW2 
							stemheight[k]=height_W2 
							W2steps[k]=1
						if (-blh[k,1] + blh[k,2]) > -(CED1 + ((WoO_T[k]-W1)/W2) * (CED2 - CED1) + (alpha * (-blh[k,0]+blh[k,2]))) and (-blh[k,0] + blh[k,2])/max((WoO_T[k]-W1),Emin) > -Eavg and (-blh[k,0] + blh[k,2])/max((WoO_T[k]-W1),Smin) < Savg:
							WoO_T[k]=WoO_T[k]+mtpervt 
							blh[k,1]=blh[k,2]
							stemheight[k]=height_W2 + (WoO_T[k]-W1) * (height_def-height_W2)/W2 
						else:
							WoO_T[k]=0 
							W2steps[k]=0 
							stemheight[k]=0 
							rnveg[k]=0 
							EPSuc[k]=0
					elif WoO_T[k] > W1+W2:
						rnveg[k]=P0 
						WoO_T[k] = -1
						W2steps[k]=0

		## print information on current timestep 
		currentDT = datetime.datetime.now()
		print ('Summer Timestep ', str(summersteps),' completed at --> ', str(currentDT),' \n')

		## Update summer season counter summersteps=summersteps+1
		summersteps=summersteps+1

	## switch between summer and winter season
	if summersteps == tstep_season and wintersteps == 0:
		WoO_T[WoO_T>1]=0	# re-set progress of WoO in every grid cell
		W2steps[WoO_T>1]=0	# re-set number of timesteps in window 2

	## winter time step calculations
	if summersteps == tstep_season and wintersteps < tstep_season: # boundary between end summer-begin winter, vegetation will be adjusted
		is_maxvalsnd=model_dfm.get_var('is_maxvalsnd')	# [-] Integral values on flow nodes
		is_sumvalsnd=model_dfm.get_var('is_sumvalsnd') 
		is_dtint=model_dfm.get_var('is_dtint')
		wd_max=is_maxvalsnd[range(ndx),2]	#maximum water level and maximum bed shear stress during the last time step
		tau_max=is_maxvalsnd[range(ndx),0]
		rslt=integrate.odeint(vegWinter, rnveg[range(ndxi)], t,h0=.05,hmax=.1) 
		rslt=integrate.odeint(vegWinter, rnveg[range(ndxi)], t)
		for veg in range(ndxi):
			if rnveg[veg]>0.1:	#decay due to bed shear stress
				rnveg[veg]=(rslt[1,[veg]])-(200/(15.0))	#New integration of cell and slow die-off due to winter: Takes 15 tsteps for max stems to disappear: 15*8/2=60 days --> Oct-Nov
			if rnveg[veg]<=0.1 and WoO_T[veg]==-1 and rnveg[veg]!=0: #if established vegetation has vanished, WOO opens and is checked again
				rnveg[veg]=0.0
				WoO_T[veg]=0 
				EPSuc[veg]=0
			elif rnveg[veg]<=0.1 and rnveg[veg]!=0:	#decay of diffused vegetation 
				rnveg[veg]=0.0
		## print information on current timestep 
		currentDT = datetime.datetime.now()
		print ('Winter Timestep ', str(wintersteps),' completed at --> ', str(currentDT),' \n')

		## Update winter season counter wintersteps=wintersteps+1
		wintersteps=wintersteps+1

	## switch between winter and summer; reset season counters 
	elif summersteps == 43 and wintersteps == 43:
		summersteps = 0
		wintersteps =0 
		model_dfm.set_var('stemheight',stemheight)
		
		## print information on current timestep 
		print ('Year ', str(year),' completed \n') 
		year=year+1
	if (i) % (1) == 0:
		rnveg_historic[:,tt]=rnveg 
		bl_historic[:,tt]=bl 
		stemheight_historic[:,tt]=stemheight 
		tt=tt+1
	
	np.savetxt('blout.txt',bl)	#export the latest bed level and vegetation values to check the output whilst running the model
	np.savetxt('rnvegout.txt',rnveg)
	
	## feed back the updated vegetation field to the FM model 
	vegupdate=np.copy(rnveg)
	vegupdate[vegupdate<2]=0	#only vegetation with densities greater than 2 stems/square meter are fed back into the DFM model
	model_dfm.set_var('rnveg',vegupdate)
	#reset the maximum bss and inundation height value counter 
	is_sumvalsnd.fill(0.0)
	is_maxvalsnd.fill(0.0) 
	is_dtint.fill(0.0)
	model_dfm.set_var('is_sumvalsnd',is_sumvalsnd) 
	model_dfm.set_var('is_maxvalsnd',is_maxvalsnd) 
	model_dfm.set_var('is_dtint',is_dtint)

model_dimr.finalize() 
print ('model finished') ## end of script
 














