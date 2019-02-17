# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 14:23:54 2019

@author: alvareo

This object is a postprocessor of the MODFLOW and MT3D simulations
"""

import os
import numpy as np
import pandas as pd
import flopy.utils.binaryfile as bf
from math import floor
from matplotlib import colors
import matplotlib.pyplot as plt

class outputdata:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,dataspace,workspace,modelname,grid):
        """
        Construction of an instance for reading data required by MODFLOW

        Parameters:
        ----------
            workspace : str
                The folder location of the dat file (private)
            name : str
                The name of the dat file, without extension      
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__dataspace=dataspace
        self.__modelname=modelname
        self.__name='output.csv'
        self.__nx=grid['nx']
        self.__ny=grid['ny']
        self.__nlay=grid['nz']
        self.__minx=grid['ox']
        self.__miny=grid['oy']
        self.__zbot=grid['oz']
        self.__delcx=grid['dx']
        self.__delry=grid['dy']
        self.__dz=grid['dz']
        self.__Lx=self.__minx+self.__delcx*self.__nx
        self.__Ly=self.__miny+self.__delry*self.__ny
        self.__Lz=self.__zbot+self.__dz*self.__nlay
    #End of constructor -------------------------------------------------------    
   
    #Begining of method -------------------------------------------------------    
    def readoutdatacsv(self): 
        """
        Method for reading the output information for posprocessing MODFLOW and MT3D results

        Parameters
        ----------        
            null
        
        Returns:
        -------        
            times : array
                Times for printing heads or concentrations
            locs: array
                coordinates of the points where piezometric curves or breakthrough
                curves will be plotted
        """         
        
        #code begins here
        flnm=self.__dataspace+'\\'+self.__name
        dfg=pd.read_csv(flnm,sep=',')
        nt=int(dfg.iloc[0,1])
        times=np.zeros(nt,dtype=np.float32)
        print('*** Reading output information, times:'+str(nt)+' ***')
        for i in range(nt):
            times[i]=dfg.iloc[i+1,1]
            print('times requested:'+str(times[i]))
        nlc=int(dfg.iloc[nt+1,1])
        print('*** Reading output information, locations:'+str(nlc)+' ***')
        locs=np.zeros((nlc,3),dtype=np.float32)
        for i in range(nlc):
            locs[i,0]=dfg.iloc[i+nt+2,0]
            locs[i,1]=dfg.iloc[i+nt+2,1]
            locs[i,2]=dfg.iloc[i+nt+2,2]
            print('locations requested: x='+str(locs[i,0])+' y='+str(locs[i,1])+' z='+str(locs[i,2]))
            
            #Convert locations on matrix index
            locs[i,0]=floor((locs[i,0]-self.__minx)/self.__delcx)
            #In MODFLOW y ans z coordinates are inverted
            locs[i,1]=floor((self.__miny+self.__delry*self.__ny-locs[i,1])/self.__delry)
            locs[i,2]=floor((self.__zbot+self.__dz*self.__nlay-locs[i,2])/self.__dz)
            print('locations requested: x='+str(locs[i,0])+' y='+str(locs[i,1])+' z='+str(locs[i,2]))
        
        return times,locs           
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def printheads(self):
        """
            Method for reading the output information for posprocessing MODFLOW and MT3D results

        Parameters:
        ----------        
            Nr, Nc : integer
                Number of rows and columns of the multimage
        
        Returns:
        -------
            null
        """
        
        #code begins here
        times,locs=self.readoutdatacsv()        
        nt=len(times)        
        hds = bf.HeadFile(os.path.join(self.__workspace,self.__modelname+'.hds'))
        extent = (self.__delcx/2., self.__Lx - self.__delcx/2., self.__Ly - self.__delry/2., self.__delry/2.)        
        
        for i in range(nt):
            axs=plt.subplot(1,1,1,aspect='equal')
            print(str(i)+" time="+str(times[i]))
            head = hds.get_data(totim=times[i])
            auxhd = head[0,:,:]
            im=axs.imshow(auxhd,interpolation='nearest',cmap=plt.cm.jet,extent=extent)            
            axs.set(title="t="+str(int(times[i]))+" d", xlabel='East (m)', ylabel='North (m)')
            plt.colorbar(im)
            plt.show()
        
    #End of method ------------------------------------------------------------    
    
    #Begining of method -------------------------------------------------------
    def printheadarray(self, Nr, Nc):
        """
        Method for reading the output information for posprocessing MODFLOW and MT3D results

        Parameters:
        ----------        
            Nr, Nc : integer
                Number of rows and columns of the multimage
        
        Returns:
        -------
            null
        """
        
        #code begins here
        
        times,locs=self.readoutdatacsv()        
        nt=len(times)       
        hds = bf.HeadFile(os.path.join(self.__workspace,self.__modelname+'.hds'))
        extent = (self.__delcx/2., self.__Lx - self.__delcx/2., self.__Ly - self.__delry/2., self.__delry/2.)
        
        fig,axs = plt.subplots(Nr,Nc,figsize=(17,9.5))        
        fig.suptitle('Piezometric heads',fontsize=14,fontweight="bold")
        images = []
        
        k=0 
        for i in range(Nr):
            for j in range(Nc):
                head = hds.get_data(totim=times[k])
                auxhd = head[0,:,:]                
                k=k+1
                images.append(axs[i,j].imshow(auxhd,interpolation='nearest',cmap=plt.cm.jet,extent=extent))
                axs[i,j].set(title="t="+str(int(times[k-1]))+" d", xlabel='East (m)', ylabel='North (m)')
                
                if(k==nt):                    
                    break                
                    
        # Find the min and max of all colors for use in setting the color scale.
        vmin = min(image.get_array().min() for image in images)
        vmax = max(image.get_array().max() for image in images)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        for im in images:
            im.set_norm(norm)
            
        fig.colorbar(images[0], ax=axs, orientation='vertical', fraction=0.02, pad=0.01)

    #Ending of method ---------------------------------------------------------
    
    
    
    