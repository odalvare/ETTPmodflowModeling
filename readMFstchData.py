# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:20:28 2018

@author: The Bodrio
"""

import os, flopy
from math import floor
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

class pospMFstch:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,workspace,name,dx,dy,dz,Lx,Ly,Lz):
        """Construction of an instance of an object for reading and displaying
        MODFLOW products after an stochastic simulation

        Parameters
        ----------
        workspace : str
            The folder location of the Modflow output (private)
        modelname : str
            The name of the Modflow model        
        
        Returns
        -------
        null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__modelname=name
        self.__delcx=dx
        self.__delcy=dy
        self.__delz=dz
        self.__Lx=Lx
        self.__Ly=Ly
        self.__Lz=Lz
        self.__nrow = floor(self.__Ly/self.__delcx)
        self.__ncol = floor(self.__Lx/self.__delry)
        self.__nlay = floor(self.__Lz/self.__delcx)
    
    #Ending of constructor ----------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def readHrealization(self,rea,t,lydisp):
        """Construction of an instance of an object for reading and displaying
        MODFLOW products after an stochastic simulation

        Parameters
        ----------
        workspace : str
            The folder location of the Modflow output (private)
        modelname : str
            The name of the Modflow model
        t : float 32
        
        Returns
        -------
        null
            
        """
        
        hds = bf.HeadFile(os.path.join(self.__workspace,self.__modelname+'.hds'))
        head = hds.get_data(totim=t)
        axuhd = head[lydisp,:,:]
        
        
    #Ending of method ---------------------------------------------------------    