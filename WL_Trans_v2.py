# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:19:43 2018

@author: The Bodrio

Stochastic model for a synthethc aquifer built upon the Walker Lake data base.
It includes flow and transport modeling on transient state for improving
control of flow and concentration variables
"""

import os, flopy
import numpy as np
#from math import floor
from Sgemsimp import Sgemsimp
from rundata import rundata
from configBAS import bldBAS
from configWELL import bldWELL
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

"""
*******************************************************************************
Assign wokspace, name of model and creation of Modflow model
"""
#Assign name and create modflow model object using work speace for storing products
directory = os.path.dirname(os.getcwd())
workspace = directory+"\\Model\\"
dataspace = directory+"\\Data\\"
figurespace = directory+"\\Figures\\"
exefile = directory+"\\Exe\\mf2005.exe"
exemt3dms = directory+'\\Exe\\mt3dms5b.exe'
fltfile = workspace+'\\mt3d_link.ftl'

#Acquisition of data for running MODFLOW
rd=rundata(dataspace)
modelname,Lx,Ly,Lz,zbot,nlay,delcx,delry,nper,steps,Lt,grid,wlay=rd.readdata()

#This is the Modflow model
mf = flopy.modflow.Modflow(modelname, exe_name=exefile, model_ws=workspace)

"""
*******************************************************************************
Discretization of the model. One layer, confined aquifer. Transient, 20 periods
each stress period has 5 intervals
"""
nrow=rd.getNrows()
ncol=rd.getNcols()
dis,nstp,perlen=rd.createdis(mf)

"""
*******************************************************************************
Next we can create a flopy object that represents the MODFLOW Basic Package.
"""
basb=bldBAS(dataspace,grid)
bas=basb.createbas(mf)

"""
*******************************************************************************
Recharge package. Uniform recharge of R=0.000001 m/d at the top layer.
"""
r = 1.0*np.ones((nrow, ncol), dtype=np.float32)/1000000.0
rcg = flopy.modflow.mfrch.ModflowRch(mf,rech=r)

"""
*******************************************************************************
Transient Well Package
"""
wellbld=bldWELL(dataspace,grid,rd)
wel=wellbld.createwel(mf)

"""
*******************************************************************************
Layer-Property Flow Package
"""
# Add LPF package to the MODFLOW model
#ipackb:A flag that is used to determine if cell-by-cell budget data should be saved
layt=[0]        #Type of layer
sg=Sgemsimp(dataspace,grid)    #101 with 50 #305:reference
hkd = 200.0*np.ones((nlay, nrow, ncol), dtype=np.float32)
                                                        #Coment the following lines when homogeneous filed is needed
hkd=sg.readsgems1rea(33)
hkd[hkd<0.001]=0.001
vkd = hkd
lpf = flopy.modflow.ModflowLpf(mf, laytyp=layt, hk=hkd, vka=vkd, ipakcb=53)

"""
*******************************************************************************
Preconditioned Conjugate Gradient Package.
"""
# PCG package for matrix computation
pcg = flopy.modflow.ModflowPcg(mf, mxiter=200, iter1=80)

"""
*******************************************************************************
Output Control
The stress period dictionary is used to set what output is saved for the corres-
ponding stress period and time step. In this case, the tuple (0, 0) means that 
stress period 1 and time step 1 for MODFLOW will have output saved. Head and 
budgets will be printed and head and budget information will be saved.
"""
# Add OC package to the MODFLOW model
#spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}
#stressperioddata={(0,2):['save head'],(1,0):['save head'],(2,1):['save head'],(3,0):['save head']}
spdata = {}
for i in range(0,nper,1):
    for j in range(0,nstp[i],1):
        spdata[i,j]=[('save head')]
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spdata, compact=True)

"""
*******************************************************************************
The Link-MT3DMS package is used to record flow information for use by MT3DMS
model : model object, the model object (of type flopy.modflow.mf.Modflow) to
which this package will be added. output_file_name : string, filename for output file
"""
# linkage to mt3dms LMT package (ATTENTION!!!)
lmt = flopy.modflow.ModflowLmt(mf, output_file_name='mt3d_link.ftl')

"""
*******************************************************************************
Writing the MODFLOW Data Files.
"""
# Write the MODFLOW model input files
mf.write_input()

"""
*******************************************************************************
Running the Modeling. 
"""
# Run the MODFLOW model
success, buff = mf.run_model()  #Flow processes are finished

"""
*******************************************************************************
Assign wokspace, name of model and creation of MT3DMS model
"""
# create mt3dms model object
mt = flopy.mt3d.Mt3dms(modflowmodel=mf, modelname=modelname, 
                       exe_name=exemt3dms , ftlfilename=fltfile,
                       model_ws=workspace)

"""
*******************************************************************************
Setting the basic parameters of transport modeling:
    prsity porosity array(nlay,nrow,ncol).
    icbund = 0 inactive, <0 constant C, >0 active, array(nlay,nr,nc).
    sconc start concentration array(nlay,nrow,ncol).
    cinact value for inactive concentration cell.
    nprs > 0 simulation saved as specified in “timprs” parameter = 0, only 
    saved at the end of simulation, < 0 saved whenever the number of transport
    steps is an even multiple of nprs
"""
# basic transport package
btn = flopy.mt3d.Mt3dBtn(mt, prsity=0.3, icbund = 1, sconc=0.0, ncomp=1, 
                         perlen = perlen, nper=nper, nstp = nstp, tsmult = 1.0, 
                         nprs = -1, nprobs = 10, cinact = -1, chkmas=True)

"""
*******************************************************************************
Objects for dealing with advection-dispersion processes
"""
# advaction package
adv = flopy.mt3d.Mt3dAdv(mt, mixelm=-1, percel=0.75)
# dispersion package
dsp = flopy.mt3d.Mt3dDsp(mt, al=1000, trpt=0.1, trpv=0.1, dmcoef=1e-09)

"""
*******************************************************************************
Setting the input of contaminants via an injection well
"""
# source/sink package
ssm_data = {}
itype = flopy.mt3d.Mt3dSsm.itype_dict()
print(itype)
print(flopy.mt3d.Mt3dSsm.get_default_dtype())

ssm_data[0] = [(0, wrow, wcol, 0.0, itype['WEL'])]
print(ssm_data)
for i in range(1,nper,1):
    if i==4:
        #wel_sp[0].append((0, wrow, wcol, pumping_rate)) # lay, row, col index, pumping rate
        #append is used for adding information on each row of the library
        ssm_data[i]=[(0, wrow, wcol, 2000.0, itype['WEL'])]    
    else:
        ssm_data[i]=[(0, wrow, wcol, 0.0, itype['WEL'])]
print(ssm_data)
ssm = flopy.mt3d.Mt3dSsm(mt, stress_period_data=ssm_data)

# matrix solver package
gcg = flopy.mt3d.Mt3dGcg(mt, cclose=1e-6)

"""
*******************************************************************************
Writing the MT3d Data Files.
"""
# write mt3dms input
mt.write_input()

"""
*******************************************************************************
Running MT3D
"""
# run mt3dms
mt.run_model()

"""
*******************************************************************************
Output of the model
"""

# plot conc
ucnobj = bf.UcnFile(os.path.join(workspace,'MT3D001.UCN'))
print(ucnobj.list_records()) # get values

wpt = ((wcol+0.5)*delry, Lx - ((wrow + 0.5)*delcx)) # origin at low upper.

times = ucnobj.get_times() # simulation time
times1 = times[round(len(times)/6.)] # 1/6 simulation time
times2 = times[round(len(times)/3.)] # 1/3 simulation time
times3 = times[round(len(times)/2.)] # 1/2 simulation time
times4 = times[round(2.*len(times)/3.)] # 2/3 simulation time
times5 = times[round(5.*len(times)/6.)] # 5/6 simulation time
times6 = times[-1] # the last simulation time

conc1 = ucnobj.get_data(totim=times1)
conc2 = ucnobj.get_data(totim=times2)
conc3 = ucnobj.get_data(totim=times3)
conc4 = ucnobj.get_data(totim=times4)
conc5 = ucnobj.get_data(totim=times5)
conc6 = ucnobj.get_data(totim=times6)

#Realization of the hydraulic conductivities from Walker Lake data base
extent = (delcx/2., Lx - delcx/2., Ly - delry/2., delry/2.)
plt.imshow(np.log10(hkd[0,:,:]), interpolation='nearest', cmap=plt.cm.jet, extent=extent, vmin=-3., vmax=5.0)
plt.colorbar()
plt.show()

#for Layer, hydraulic heads
totaltimes = np.cumsum(perlen)
print(totaltimes)
plt.subplot(1,1,1,aspect='equal')
levels = np.arange(0.,12.,0.5)
hds = bf.HeadFile(os.path.join(workspace,modelname+'.hds'))

for i in range(0,totaltimes.size,1):
    head = hds.get_data(totim=totaltimes[i])
    axuhd = head[0,:,:]
    plt.imshow(axuhd, interpolation='nearest', cmap=plt.cm.jet, extent=extent)
    plt.colorbar()
    plt.show()
    
# conc 100 day
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc1) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times1)
plt.show()

# conc in 500 days
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc2) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times2)
plt.show()

# conc in 1000 days
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc3) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times3)
plt.show()

# conc 100 day
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc4) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times4)
plt.show()

# conc in 500 days
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc5) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times5)
plt.show()

# conc in 1000 days
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
#lc = modelmap.plot_grid() # grid
cs = modelmap.plot_array(conc6) # head contour
plt.colorbar(cs) # colrobar
plt.plot(wpt[0],wpt[1],'ro')
plt.title('C  %g day' % times6)
plt.show()

