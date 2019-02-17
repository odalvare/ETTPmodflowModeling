# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:12:47 2019

@author: alvareo

This clas assemble the wells dictionary required for the calling of
MODFLOW WELL module constructor

"""

import pandas as pd
import numpy as np
import flopy
from math import floor
from rundata import rundata

class bldWELL:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,workspace,gridin,df,csim):
        """
        Construction of an instance for reading well data required by MODFLOW

        Parameters
        ----------
            workspace : string
                The folder location of the sgems output (private)
            grid : dictionary
                Contains the information of the domain's geometry
            df : rundata object
                Object containing the running information of the MODFLOW model
            csim : boolean
                Flag for declaring concentration simulation via MT3D. If True
                the program reads the dictionaries with concentration injection
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__grid=gridin
        self.__dflst=pd.read_csv(workspace+'WELL\\'+'filelist.csv',sep=',')
        print('*--- Preparing simulation, reading csv wells ---*')
        self.__dfwells=pd.read_csv(workspace+'WELL\\'+'arraydata.csv',sep=',')
        self.__nwells=self.__dflst.shape[0]
        self.__df=df
        self.__dfclst={}
        self.__dfcwells={}
        self.__ncwells=0
        if(csim==True):
            self.buildCwells()        
        
    #End of constructor -------------------------------------------------------
    
     #Begining of constructor --------------------------------------------------
    def buildCwells(self):
        """
        Read he files required for WELL package of the MT3D code for transport
        moeling with MODFLOW

        Parameters
        ----------
            null
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__dfclst=pd.read_csv(self.__workspace+'WELL\\'+'fileclist.csv',sep=',')
        print('*--- Preparing transport simulation, reading csv concentration wells ---*')
        self.__dfcwells=pd.read_csv(self.__workspace+'WELL\\'+'arraycdata.csv',sep=',')
        self.__ncwells=self.__dflst.shape[0]
        
    #End of constructor -------------------------------------------------------
    
    #Beginning of method ------------------------------------------------------
    def getDictWells(self):
        """
        This method constructs the dictionary of wells for the MODFLOW simulation
        Each compoent of the dicxtionary represents the forcings acting on
        each forcing period, thus len(dict)=nper (on rundata object).

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            wel_sp : dictionary of arrays
                3dictionary of wells for the MODFLOW simulation
                
        .. note:
            * in folder dataspace\\WELLS\\ must be located files named:
                filelist.csv, lis of wells and their characteristics
                arraydata.csv, pumping data for each well (with date and as many
                rowws as forcint periods of the model)
        """
        #Method begins here
        #nx=self.__grid['nx']   #From the geometry in grid
        ny=self.__grid['ny']
        nz=self.__grid['nz']
        minx=self.__grid['ox']
        miny=self.__grid['oy']
        minz=self.__grid['oz']
        rx=self.__grid['dx']
        ry=self.__grid['dy']
        rz=self.__grid['dz']
        
        # well package
        # Remember to use zero-based layer, row, column indices!
        lcoordw=np.zeros((self.__nwells,3),dtype=np.int32)
        for i in range (self.__nwells):
            lcoordw[i,0]=floor((self.__dflst.iloc[i,3]-minx)/rx)
            #In MODFLOW y ans z coordinates are inverted
            lcoordw[i,1]=floor((miny+ry*ny-self.__dflst.iloc[i,4])/ry)
            lcoordw[i,2]=floor((minz+rz*nz-self.__dflst.iloc[i,5])/rz)
        
        nper=self.__df.getForcPer()
        wel_sp = {}        
        for i in range(nper):
            lst=[]
            for j in range(self.__nwells):
                pumping_rate=self.__dfwells.iloc[i+1,j+1]
                lst.append( [lcoordw[j,2], lcoordw[j,1], lcoordw[j,0], pumping_rate] )
            wel_sp[i]=lst
        print(wel_sp)
        
        print('*--- Succesfull reading of wells ---*')
        
        return wel_sp
    #End of method ------------------------------------------------------------
    
    #Beginning of method ------------------------------------------------------
    def getDictCWells(self,itype):
        """
        This method constructs the dictionary of wells for the MODFLOW simulation
        Each compoent of the dicxtionary represents the forcings acting on
        each forcing period, thus len(dict)=nper (on rundata object).

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            wel_sp : dictionary of arrays
                3dictionary of wells for the MODFLOW simulation
                
        .. note:
            * in folder dataspace\\WELLS\\ must be located files named:
                filelist.csv, lis of wells and their characteristics
                arraydata.csv, pumping data for each well (with date and as many
                rowws as forcint periods of the model)
        """
        #Method begins here
        #nx=self.__grid['nx']   #From the geometry in grid
        ny=self.__grid['ny']
        nz=self.__grid['nz']
        minx=self.__grid['ox']
        miny=self.__grid['oy']
        minz=self.__grid['oz']
        rx=self.__grid['dx']
        ry=self.__grid['dy']
        rz=self.__grid['dz']
        
        # well package
        # Remember to use zero-based layer, row, column indices!
        lcoordw=np.zeros((self.__ncwells,3),dtype=np.int32)
        for i in range (self.__ncwells):
            lcoordw[i,0]=floor((self.__dfclst.iloc[i,3]-minx)/rx)
            #In MODFLOW y ans z coordinates are inverted
            lcoordw[i,1]=floor((miny+ry*ny-self.__dfclst.iloc[i,4])/ry)
            lcoordw[i,2]=floor((minz+rz*nz-self.__dfclst.iloc[i,5])/rz)
        
        nper=self.__df.getForcPer()
        ssm_data = {}
        print('Number of conc periods='+str(nper))        
        for i in range(nper):
            lst=[]
            for j in range(self.__ncwells):
                conc_rate=self.__dfcwells.iloc[i+1,j+1]
                lst.append( [ lcoordw[j,2], lcoordw[j,1], lcoordw[j,0], conc_rate, itype['WEL'] ] )
            ssm_data[i]=lst
        print(ssm_data)
        
        print('*--- Succesfull reading of concentration wells ---*')
        
        return ssm_data
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getDataDict(self): 
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
        return self.__dflst,self.__dfwells
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNoWells(self): 
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
        return self.__nwells
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def createwel(self,mf):
        """
        Build *.well file for a MODFLOW simulation

        Parameters:
        ----------        
            mf: modflow objecto from flopy
                Object to atttach the bas file
        
        Returns:
        -------
            wel : flopy.modflow.ModflowWel object
                wel object attached to the mf object of a MODFLOW simulation
        """
        wel_sp=self.getDictWells()
        wel = flopy.modflow.ModflowWel(mf, stress_period_data=wel_sp)
        return wel
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def createssm(self,mt,itype):
        """
        Build *.ssm file for a MT3D simulation

        Parameters:
        ----------        
            mt: mt3d object from flopy
                Object to atttach the bas file
        
        Returns:
        -------
            ssm : flopy.mt3d.Mt3dSsm
                concentration ssm object attached to the MT3D simulation
        """
        ssm_data=self.getDictCWells(itype)
        ssm = flopy.mt3d.Mt3dSsm(mt, stress_period_data=ssm_data)
        return ssm
    #End of method ------------------------------------------------------------
    
    