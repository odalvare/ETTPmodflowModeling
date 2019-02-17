# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 09:04:05 2018

@author: The Bodrio
"""

import numpy as np
import pandas as pd

class Sgemsimp:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,workspace,grid):
        """Construction of an instance of an object for reading sgmes output files

        Parameters
        ----------
            workspace : string
                The folder location of the sgems output (private)
            grid : dictionary
                Contains the information of the domain's geometry
        
        Returns
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__grid=grid
        
    #Ending of constructor ----------------------------------------------------
        
    #Begining of method -------------------------------------------------------
    def readsgems1rea(self,rea):
        """
        Method for reading the sgems csv file for one realization produced
        by a simulation (or estimation) algorithm in sgems

        Parameters:
        ----------        
            rea : integer
                Realization used for the transient deterministic simulation
        
        Returns:
        -------
            datar : float array
                3D array containing the data of the field in order (col->row->lay)
        """
        #Method begins here, use numpy np.getfromtxt for extracting the column of realization
        print('*--- Preparing simulation, reading K info ---*')
        flnm=self.__workspace+'sims.csv'
        print('File with K data='+flnm)
        cols=['real'+str(rea)]
        reasercol = pd.read_csv(flnm,sep=',',usecols=cols,squeeze=True)
        count_row = reasercol.shape[0]                                         # number of row count, number of simulated data
        print('K values='+str(count_row))
        readata = np.reshape(reasercol, (count_row), order='F')                #Delete zero entrys and replace them with min
        readata[readata==0.]=np.nan
        minv=np.amin(readata)
        readata[np.isnan(readata)]=minv 
        
        cols=self.__grid['nx']
        rows=self.__grid['ny']
        lays=self.__grid['nz']
        datar=np.zeros((lays,rows,cols),dtype='float32')                       #Init array of data to use in MODFLOW
        
        cont=0
        for k in range(0,lays,1):                                              #Reading of data
            for i in range(0,rows,1):
                for j in range(0,cols,1):
                    datar[lays-k-1,rows-i-1,j]=readata[cont]                  #Invert rows and layers, origin is on southrn-western corner
                    cont=cont+1
        
        return datar
    
    #End of method ------------------------------------------------------------