# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 10:54:57 2019

@author: alvareo

This clas assemble the ibopund and strt arrays required for the calling of
MODFLOW BAS module constructor

"""

from ascii import ascii
import numpy as np
import flopy

class bldBAS:
    
    #Beginning of constructor -------------------------------------------------
    def __init__(self,workspace,gridin):
        """
        Construction of an instance for reading data required by MODFLOW

        Parameters
        ----------
            workspace : string
                The folder location of the sgems output (private)
            grid : dictionary
                Contains the information of the domain's geometry      
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__grid=gridin
        
    #End of constructor -------------------------------------------------------
    
    #Beginning of method ------------------------------------------------------
    def getibounds(self):
        """
        Method for assamble the ibound integer array cell with values 
        (-1:constant head, 0:inactive cells, 1:activecells) from ascii raster 
        text files

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            ibound : int32 array
                3D array containing the ibound of the field (lay,row,col)
                
        .. note:
            * in folder dataspace\\IBOUND\\ must be as many ascii files as
            layers of the model. Otherwise, the execution crashes.
        """
        #Method begins here
        ncol=self.__grid['nx']              #From the geometry in grid
        nrow=self.__grid['ny']
        nlay=self.__grid['nz']
        minx=self.__grid['ox']
        miny=self.__grid['oy']
        rx=self.__grid['dx']
        ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)    #layers,rows,cols
        
        print('*--- Preparing simulation, reading ibound ---*')
        
        for i in range(nlay):
            print('Layer '+str(i+1))
            asc=ascii(self.__workspace+'IBOUND\\',"ibound"+str(i+1))
            colP,filP,xllP,yllP,dxP,ndataP=asc.header()
            asc.checkprmtrs(nc=ncol,nf=nrow,xi=minx,yi=miny,r=rx,miss=-9999.0)
            mtrxib=asc.readmtrx(2)
            ibound[i,:,:]=mtrxib[:,:]
        
        print('*--- Succesfull reading of ibound ---*')
        
        return ibound
    #End of method ------------------------------------------------------------
    
    #Beginning of method ------------------------------------------------------
    def getstrt(self):
        """
        Method for assamble the strt float array (imposed and initial heads)
        from ascii raster  text files

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            strt : float32 array
                3D array containing the strt of the field (lay,row,col)
                
        .. note:
            * in folder dataspace\\STRT\\ must be as many ascii files as
            layers of the model. Otherwise, the execution crashes.
        """
        #Method begins here
        ncol=self.__grid['nx']              #From the geometry in grid
        nrow=self.__grid['ny']
        nlay=self.__grid['nz']
        minx=self.__grid['ox']
        miny=self.__grid['oy']
        rx=self.__grid['dx']
        strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)    #layers,rows,cols
        
        print('*--- Preparing simulation, reading ibound ---*')
        
        for i in range(nlay):
            print('Layer '+str(i+1))
            asc=ascii(self.__workspace+'STRT\\',"strt"+str(i+1))
            colP,filP,xllP,yllP,dxP,ndataP=asc.header()
            asc.checkprmtrs(nc=ncol,nf=nrow,xi=minx,yi=miny,r=rx,miss=-9999.0)
            mtrxib=asc.readmtrx(2)
            strt[i,:,:]=mtrxib[:,:]
        
        print('*--- Succesfull reading of strt ---*')
        
        return strt
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def createbas(self,mf):
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
        ibound=self.getibounds()
        strt=self.getstrt()
        bas=flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
        return bas
    #End of method ------------------------------------------------------------
        