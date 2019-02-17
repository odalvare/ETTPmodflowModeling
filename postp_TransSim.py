# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:39:53 2019

@author: The Bodrio

Extraction and post processing of transport modeling using MT3D via flopy
"""

import os, flopy
import numpy as np
from math import floor
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

"""
*******************************************************************************
Assign wokspace, name of model and creation of Modflow model
"""
#Assign name and create modflow model object using work speace for storing products
directory = os.path.dirname(os.getcwd())
modelname = 'WalkerLake'
workspace = directory+"\Model"
#workspace1 = directory + '\Backup\NewGradient'
dataspace = directory+"\Data"
figurespace = directory+"\Figures"
exefile = directory+"\Exe\mf2005.exe"
exemt3dms = directory+'\Exe\mt3dms5b.exe'
fltfile = workspace+'\mt3d_link.ftl'
nmconcfile = 'MT3D001.UCN'

"""
*******************************************************************************
Output of the model
"""

# Model domain and grid definition
prsity = 0.3
Lx = 220.
Ly = 140.
ztop = 0.
zbot = -50.
nlay = 1
delcx = 1.
delry = 1.
nrow = floor(Ly/delcx)
ncol = floor(Lx/delry)
delv = (ztop - zbot) / nlay
widhtly = [50.]
botm = np.linspace(ztop, zbot, nlay + 1)
for i in range (1,len(botm),1):
    botm[i]=botm[i-1]-widhtly[i-1]

extent = (delcx/2., Lx - delcx/2., Ly - delry/2., delry/2.)

# Extract concentration for each simulation time
ucnobj = bf.UcnFile(os.path.join(workspace,nmconcfile))
print(ucnobj.list_records())    # get values
times = ucnobj.get_times()      # simulation time
print("Simulation times="+str(len(times)))

#plt.subplot(1,1,1,aspect='equal')
#levels = np.arange(0.,100.,20)

intconc = np.zeros( (2,len(times)), dtype=np.float32 ) #Array for storing the integral of concentrations
intmass = np.zeros( (2,len(times)), dtype=np.float32 ) #Array for storing the accumulated mass
btc1 = np.zeros( (2,len(times)), dtype=np.float32 ) #Array for storing the breakthrough curve 1
btc2 = np.zeros( (2,len(times)), dtype=np.float32 ) #Array for storing the breakthrough curve 1
btc3 = np.zeros( (2,len(times)), dtype=np.float32 ) #Array for storing the breakthrough curve 1
refmlc = np.zeros( (2,2), dtype=np.float32 )

for t in range(0,len(times),1):
    time=times[t]
    print("Extracting concentrations in time="+str(time))
    conc = ucnobj.get_data(totim=time) #Extract a 3D array containing concentration in time t
    aux = conc[0,:,:]                #Puts concentration in a 2D array    
    
    #Calculate average concentration on each time for the aquifer
    avc=0.
    am=0.
    cont=0
    for k in range(0,nlay,1):
        for i in range(0,nrow,1):
            for j in range(0,ncol,1):
                if aux[i,j]<0.0001:
                    aux[i,j]=0.0001
                if ((i==100)and(j==80)):
                    btc1[0,t]=time
                    btc1[1,t]=aux[i,j]
                if ((i==90)and(j==150)):
                    btc2[0,t]=time
                    btc2[1,t]=aux[i,j]
                if ((i==80)and(j==210)):
                    btc3[0,t]=time
                    btc3[1,t]=aux[i,j]
                avc=avc+aux[i,j]*delcx*delry*delv*prsity*(1000.0/len(times))
                #am=am+aux[i,j]*delcx*delry*delv*prsity*(1000.0/len(times))
                cont=cont+1                
    
    avc=(avc)/1000000.
    #am=am/cont
    intconc[0,t]=time
    intconc[1,t]=avc
    intmass[0,t]=time
    intmass[1,t]=am
    print("Average concentration="+str(avc)+"  Accumulated mass="+str(avc)) 

refmlc[0,0]=0.
refmlc[0,1]=1000.
refmlc[1,0]=-1.
refmlc[1,1]=-1.
    
fig, ax = plt.subplots()
ax.plot(intconc[0,:], intconc[1,:])
ax.set(xlabel='time (day)', ylabel='Mass (Mg)',title='Accumulated mass graph')
ax.grid()

figa, axa = plt.subplots()
axa.plot(btc1[0,:], np.log10(btc1[1,:]))
axa.plot(refmlc[0,:],refmlc[1,:])
axa.set(xlabel='time (day)', ylabel='log concentration (g/cub.m)',title='Breakthrough curve i=100 j=80')#, vmin=-4.0, vmax=2.0)
axa.grid()

figb, axb = plt.subplots()
axb.plot(btc2[0,:], np.log10(btc2[1,:]))
axb.plot(refmlc[0,:],refmlc[1,:])
axb.set(xlabel='time (day)', ylabel='log concentration (g/cub.m)',title='Breakthrough curve i=90 j=150')#, vmin=-4.0, vmax=2.0)
axb.grid()

figc, axc = plt.subplots()
axc.plot(btc3[0,:], np.log10(btc3[1,:]))
axc.plot(refmlc[0,:],refmlc[1,:])
axc.set(xlabel='time (day)', ylabel='log concentration (g/cub.m)',title='Breakthrough curve i=80 j=210')#, vmin=-4.0, vmax=2.0)
axc.grid()

#dispt=1000.
#conc = ucnobj.get_data(totim=dispt) #Extract a 3D array containing concentration in time t
#aux = np.log10(conc[0,:,:])                #Puts concentration in a 2D array
#aux2 = np.zeros((nrow, ncol), dtype=np.float32)  
#for i in range(0,nrow,1):
#    for j in range(0,ncol,1):
#        if ( (aux[i,j]>3.) or (aux[i,j]<-3.) ):
#            aux2[i,j]=np.nan
#        else:
#            aux2[i,j]=aux[i,j]            
#
#levels = [-1]        
#plt.imshow(aux2, interpolation='nearest', cmap=plt.cm.jet, extent=extent, vmin=-3., vmax=3.0)
#plt.title("Log concentration at time "+str(dispt))
#plt.colorbar()
#plt.contour(aux2,levels)
#plt.show()