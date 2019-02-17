# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 08:33:14 2019

@author: alvareo

For testing objects and executions

"""

import os, flopy
from runMF import runMF
from runMT3D import runMT3D
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
#from outputdata import outputdata
#from configWELL import bldWELL
#from rundata import rundata

directory = os.path.dirname(os.getcwd())
workspace = directory+"\\Model\\"
dataspace = directory+"\\Data\\"
outspace = directory+'\\Output\\'
figurespace = directory+"\\Figures\\"
exefile = directory+"\\Exe\\mf2005.exe"
exemt3dms = directory+'\\Exe\\mt3dms5b.exe'
fltfile = outspace+'mt3d_link.ftl'

rea=38
thrs=0.001

mfexe=runMF(exefile,workspace,dataspace,rea,thrs,csim=True)
mfexe.runMODFLOW(fltfile)
grid=mfexe.getGrid()
modelname=mfexe.getMfName()
mf=mfexe.getModel()

mtexe=runMT3D(mf,mfexe,exemt3dms,modelname,workspace,dataspace,fltfile)
mtexe.runMT3D()

nx=grid['nx']
ny=grid['ny']
nlay=grid['nz']
minx=grid['ox']
miny=grid['oy']
zbot=grid['oz']
delcx=grid['dx']
delry=grid['dy']
dz=grid['dz']
Lx=minx+delcx*nx
Ly=miny+delry*ny
Lz=zbot+dz*nlay
wcol=25
wrow=70

## plot conc
ucnobj = bf.UcnFile(os.path.join(workspace,'MT3D001.UCN'))
print(ucnobj.list_records()) # get values
#
wpt = ((wcol+0.5)*delry, Lx - ((wrow + 0.5)*delcx)) # origin at low upper.
#
times = ucnobj.get_times() # simulation time
times1 = times[round(len(times)/6.)] # 1/6 simulation time
times2 = times[round(len(times)/3.)] # 1/3 simulation time
times3 = times[round(len(times)/2.)] # 1/2 simulation time
times4 = times[round(2.*len(times)/3.)] # 2/3 simulation time
times5 = times[round(5.*len(times)/6.)] # 5/6 simulation time
times6 = times[-1] # the last simulation time
#
conc1 = ucnobj.get_data(totim=times1)
conc2 = ucnobj.get_data(totim=times2)
conc3 = ucnobj.get_data(totim=times3)
conc4 = ucnobj.get_data(totim=times4)
conc5 = ucnobj.get_data(totim=times5)
conc6 = ucnobj.get_data(totim=times6)

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
