# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 10:58:47 2019

@author: alvareo

Object for extracting the main data required by MODFLOW execution and display
of results. Also configures the discretization of the model in MODFLOW and M3T3D

"""

import numpy as np
from math import floor
import flopy

class rundata:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,workspace):
        """
        Construction of an instance for reading data required by MODFLOW

        Parameters:
        ----------
            workspace : string
                The folder location of the dat file (private)
            name : string
                The name of the dat file, without extension      
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__name='runMF'
        self.__nmodel=""
        self.__minx=float(0.0)
        self.__Lx=float(0.0)
        self.__miny=float(0.0)
        self.__Ly=float(0.0)
        self.__Lz=float(0.0)
        self.__Lt=float(0.0)
        self.__zbot=float(0.0)
        self.__nlay=int(1)
        self.__nx=int(1)
        self.__ny=int(1)
        self.__delcx=int(1)
        self.__delry=int(1)
        self.__delz=int(1)
        self.__nper=int(20)
        self.__nsteps=int(1)
        self.__grid={}
        self.__wlay=[]
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------    
    def readdata(self,verbose=False): 
        """
        Method for reading the input data for a MODFLOW simulation

        Parameters
        ----------        
            null
        
        Returns
        -------        
            (nmodel,Lx,Ly,Lz,zbot,nlay,delcx,delry,npere,steps) : tuple
                required data for running MODFLOW.           
        """
        #code begins here
        
        nmfl=self.__workspace+self.__name+'.dat'
        
        try:
            
             with open(nmfl, 'rb') as in_file:
                # A dictionary to contain info about the grid                
                if verbose :
                    print(('\n    Reading file: "{0}"'.format(nmfl)))
                
                print('*--- Preparing simulation, reading grid info ---*')
                line=in_file.readline()
                line_split=line.split()
                self.__nmodel=str(line_split[1])
                self.__nmodel=self.__nmodel[2:-1]
                print('model name='+self.__nmodel)                             #Model name
                
                line=in_file.readline()
                line_split=line.split()
                self.__minx=float(line_split[1])
                self.__grid['ox']=self.__minx
                print('Xorigin='+str(self.__minx))                             #Model Xorigin
                
                line=in_file.readline()
                line_split=line.split()
                self.__Lx=float(line_split[1])
                print('Xlength='+str(self.__Lx))                               #Model Lx
                
                line=in_file.readline()
                line_split=line.split()
                self.__miny=float(line_split[1])
                self.__grid['oy']=self.__miny
                print('Yorigin='+str(self.__miny))                             #Model Yorigin
                
                line=in_file.readline()
                line_split=line.split()
                self.__Ly=float(line_split[1])
                print('Ylength='+str(self.__Ly))                               #Model Ly
                
                line=in_file.readline()
                line_split=line.split()
                self.__zbot=float(line_split[1])
                self.__grid['oz']=self.__zbot
                print('Zorigin='+str(self.__zbot))                             #Model Zorigin
                
                line=in_file.readline()
                line_split=line.split()
                self.__Lz=float(line_split[1])
                print('Zlength='+str(self.__Lz))                               #Model Lz
                
                line=in_file.readline()
                line_split=line.split()
                self.__nlay=int(line_split[1])
                self.__grid['nz']=self.__nlay
                print('Nlayers='+str(self.__nlay))                             #Model Layers
                
                line=in_file.readline()
                line_split=line.split()
                self.__delcx=float(line_split[1])
                self.__grid['dx']=self.__delcx
                print('Deltax='+str(self.__delcx))                             #Model DeltaX
                
                line=in_file.readline()
                line_split=line.split()
                self.__delry=float(line_split[1])
                self.__grid['dy']=self.__delry
                print('Deltay='+str(self.__delry))                             #Model DeltaY
                
                self.__nx=floor(self.__Lx/self.__delcx)
                self.__grid['nx']=self.__nx
                print('#Cols='+str(self.__nx))                                 #Model Number of columns
                
                self.__ny=floor(self.__Ly/self.__delry)
                self.__grid['ny']=self.__ny
                print('#Rows='+str(self.__ny))                                 #Model Number of rows
                
                self.__delz=floor(self.__Lz/self.__nlay)
                self.__grid['dz']=self.__delz
                print('DeltaZ='+str(self.__delz))                              #Model DeltaX
                
                line=in_file.readline()
                line_split=line.split()
                self.__nper=int(line_split[1])                                 #Simulation periods
                print('#Pers='+str(self.__nper))
                
                line=in_file.readline()
                line_split=line.split()
                self.__nsteps=int(line_split[1])                               #Simulation forcing steps
                print('#Steps='+str(self.__nsteps))
                      
                line=in_file.readline()
                line_split=line.split()
                self.__Lt=float(line_split[1])                                 #Simulation horizon in days
                print('#Steps='+str(self.__Lt))
                
                #Read the width of the layers
                line=in_file.readline()
                line=in_file.readline()
                self.__wlay=np.zeros((self.__nlay),dtype=np.float32)
                    
                for i in range (self.__nlay):
                    line=in_file.readline()
                    self.__wlay[i]=float(line)
                
                self.__grid['type'] = 'points'
                self.__grid['points']  = self.__nx*self.__ny*self.__nlay
                
                in_file.close()

                print('*--- Preparing simulation, grid info succesfully read ---*')                
            
        except IOError:
            print(('    Error reading file "{0}"'.format(nmfl)))
            print('    Check if the file exists...')

        print('*---Succesfull reading of the execution data in file '+self.__name+'.dat---*')
        
        return  (self.__nmodel,self.__Lx,self.__Ly,self.__Lz,self.__zbot,
                self.__nlay,self.__delcx,self.__delry,self.__nper,self.__nsteps,
                self.__Lt,self.__grid,self.__wlay)           
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNcols(self): 
        """
        Get the number of columns
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nx: int
                Number of columns on the grid
        """
        #code begins here
        return self.__nx
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNrows(self): 
        """
        Get the number of row
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            ny: int
                Number of rows on the grid
        """
        #code begins here
        return self.__ny
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNlays(self): 
        """
        Get the number of layers
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nx: int
                Number of layes on the grid
        """
        #code begins here
        return self.__nz
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getBottoms(self): 
        """
        Get the altitudes of the bottoms of the layers
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            botm: array float
                Bottoms of the layers
        """
        #code begins here
        ztop=self.__zbot+self.__Lz
        widhtly=self.__wlay[0:self.__nlay]
        botm = np.linspace(ztop,self.__zbot,self.__nlay+1)
        for i in range (1,len(botm),1):
            botm[i]=botm[i-1]-widhtly[i-1]
        return botm
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getForcPer(self): 
        """
        Get the number of forcing periods of the simulation
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nper: int
                Number of forcing periods of the simulation
        """
        #code begins here        
        return self.__nper
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getSteps(self): 
        """
        Get the number of steps within a forcing period
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nsteps: int
                Number of steps within a forcing period
        """
        #code begins here
        return self.nsteps
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def createdis(self,mf):
        """
        Build *.bas file for a MODFLOW simulation

        Parameters:
        ----------        
            mf: modflow objecto from flopy
                Object to atttach the bas file
        
        Returns:
        -------
            bas : flopy.modflow.ModflowBas object
                bas object attached to the mf object of a MODFLOW simulation
        """
        # Model domain and grid definition
        ztop=self.__zbot+self.__Lz
        botm=self.getBottoms()
        
        perlen=np.ones(self.__nper, dtype=np.float32)
        perlen=[self.__Lt/float(self.__nper) for x in perlen]                  #Stress period duration in days
        #perlen=[x*86400 for x in perlen] #                                      en segundos (SI)
        nstp=np.ones(self.__nper, dtype=np.int32)
        for i in range(0,self.__nper,1):                               #Number of time steps in each stress period
            nstp[i]=self.__nsteps
            
        #With this information, we can now create the flopy discretization object by entering the following:
        #Default time units: lenuni=2:meters and itmuni=4:days
        dis = flopy.modflow.ModflowDis(mf, self.__nlay, self.__ny, self.__nx,
                                       delr=self.__delry, delc=self.__delcx,
                                       top=ztop, botm=botm[1:], nper=self.__nper,
                                       perlen=perlen, nstp=nstp, steady=False,
                                       itmuni=4, lenuni=2)
        return (dis,nstp,perlen)
    #End of method ------------------------------------------------------------
        
    